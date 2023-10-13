/* Linear time Closed itemset Miner for Frequent Itemset Mining problems */
/* 2004/4/10 Takeaki Uno */
/* This program is available for only academic use.
   Neither commercial use, modification, nor re-distribution is allowed */

#ifndef _lcm_c_
#define _lcm_c_

#include<time.h>
#include<stdio.h>
#include"lib_e.c"
#include"lcm_var.c"
#define LCM_PROBLEM LCM_CLOSED
#include"trsact.c"
#include"lcm_io.c"
#include"lcm_init.c"
#include"lcm_lib.c"

/* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
#include"transaction_keeping.c"
/* END OF MODIFICATIONS */

#include"fastyb.c"

/* MODIFICATIONS TO KEEP TRACK OF EXECUTION TIME AND MEMORY CONSUMPTION */
#include"time_keeping.c"
/* END OF MODIFICATIONS */

/* FUNCTION DECLARATIONS OF ORIGINAL LCM SOURCE */
void LCMclosed_BM_iter(int item, int m, int pmask);

void LCMclosed_BM_recursive (int item, int mask, int pmask){
  int i;

  if((i=LCM_Ofrq[0]) >= LCM_th){
    if((LCM_BM_pp[1]&pmask)==1 && LCM_BM_pt[1]==0)
        LCM_print_last(LCM_Op[0], i);
//     else printf ("11clo item%d %d: %d : %x & %x = %x,   %d:  prv%d pprv%d \n", LCM_Op[0], LCM_itemsett, LCM_frq, pmask, LCM_BM_pp[1], LCM_BM_pp[1] & pmask, LCM_BM_pt[1], LCM_prv, LCM_pprv);
   // pruning has to be here
  }
  LCM_BM_weight[1] = LCM_Ofrq[0] = 0;
  LCM_Ot[0] = LCM_Os[0];
  /* MODIFICATIONS TO KEEP TRACK OF TRANSACTIONS */
  BM_TRANS_LIST_EMPTY(1);
  /* END OF MODIFICATIONS */
  for(i=1; i<item; i++){
    if(LCM_Ofrq[i] >= LCM_th) LCMclosed_BM_iter(i, mask, pmask);
    else{
    	if(LCM_Ofrq[i]>0) {
    		//printf("Item %d is infrequent (%d,%d) and will not have LCMclosed_BM_iter,\n",i,LCM_Ofrq[i],LCM_th);
    		LCM_BM_occurrence_delete(i);
    		//printf("Problem fixed!\n");
    	}
    }
  }
  	 // Maybe an else if (LCM_Ofrq[i] > 0) LCM_BM_occurrence_delete(i) doesn't hurt
}



/*************************************************************************/
/* LCMclosed iteration (bitmap version ) */
/* input: T:transactions(database), item:tail(current solution) */
/*************************************************************************/
void LCMclosed_BM_iter(int item, int m, int pmask){
  int mask, it = LCM_itemsett, ttt;

  LCM_frq = LCM_Ofrq[item];
  pmask &= BITMASK_31[item];
  if((ttt = LCM_BM_closure(item, pmask)) > 0){
   // pruning has to be here
	  //printf ("BMclo %d item%d it%d frq%d,  prv%d pprv%d::  ttt=%d,%d pmask%x\n", item, LCM_Op[item], LCM_itemsett, LCM_frq, LCM_prv, LCM_pprv, ttt,LCM_Op[ttt], pmask );
    LCM_BM_occurrence_delete(item);
    return;
  }
  LCM_iters++;
  BUF_reset(&LCM_B);
  LCMclosed_BM_occurrence_deliver_(item, m);
  LCM_additem(LCM_Op[item]);
  mask = LCM_BM_rm_infreq(item, &pmask);

  bm_process_solution_ap_perm(LCM_frq,item,&mask); // PAOLO
  LCM_solution();
  /* MODIFICATION FOR FAST WY ALGORITHIM */
  if(LCM_frq != current_trans.siz){
	  printf("LCM_frq=%d, current_trans.siz=%d\n",LCM_frq,current_trans.siz);
  }

  /* END OF MODIFICATION */

  /* MODIFICATION TO KEEP TRACK OF TRANSACTIONS */
  //print_current_trans();
  BM_CURRENT_TRANS_EMPTY();
  /* END OF MODIFICATION */

  LCMclosed_BM_recursive(item, mask, pmask);
  while(LCM_itemsett>it) LCM_delitem();
  BUF_clear(&LCM_B);
}




/***************************************************************/
/* iteration of LCM ver. 2 */
/* INPUT: T:transactions(database), item:tail of the current solution */
/*************************************************************************/
// LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
int LCMclosed_iter(ARY *T, int item, int prv, TRANS_LIST *trans_list){
  ARY TT;
  int i, ii, e, ee, n, js=LCM_jump.s, it=LCM_itemsett, mask;
  int flag=-1, perm[LCM_BM_MAXITEM], pmask = 0xffffffff;
  QUEUE_INT *q;

  LCM_jump.s = LCM_jump.t;
  LCM_iters++;
  LCM_additem(item);
  LCM_frq = LCM_Ofrq_[item];
  LCM_prv = item;
  LCM_pprv = prv;

  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  TRANS_LIST mk_trans_list, shrink_trans_list;
  /* END OF MODIFICATIONS */

  //printf( " Ot %d %d (%d,%d)\n", LCM_Ot[item]-LCM_Os[item], LCM_frq, item,prv);
  n = LCM_freq_calc(T, item, LCM_Eend-1);
  if(prv >= 0) LCM_Ofrq[prv] = 0;
  LCM_Ofrq[item] = 0;
  ii = LCM_jump_rm_infreq(item);
  LCM_jumpt = LCM_jump.t;
  if(ii > item){
    flag = ii;
    //printf ("###clo item%d %d %d: %d\n", item, ii, LCM_Ofrq[ii], LCM_Ofrq_[ii]);
    goto END2;
  }  /* itemset is not closed */

  BUF_reset(&LCM_B);
  LCM_partition_prefix(item);

  if(QUEUE_LENGTH(LCM_jump)==0){
    LCM_Ofrq[item] = 0;
    LCMclosed_BM_occurrence_deliver_first(item, T, trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
    mask = LCM_BM_rm_infreq(LCM_BM_MAXITEM, &pmask);

    ary_process_solution_ap_perm(LCM_frq, trans_list, item, &mask); // PAOLO
    LCM_solution();
    /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
    //print_transaction_list(trans_list,item);
    /* END OF MODIFICATIONS */
    LCMclosed_BM_recursive(LCM_BM_MAXITEM, mask, pmask);
    BUF_clear(&LCM_B);
    goto END2;
  }

  LCM_BM_occurrence_deliver_first(item, T);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  QUEUE_FE_LOOP_(LCM_jump, i, ii) LCM_Ofrq[ii] = 0;
  mask = LCM_BM_rm_infreq(LCM_BM_MAXITEM, &pmask);

  ary_process_solution_ap_perm(LCM_frq, trans_list, item, &mask); // PAOLO
  LCM_solution();
  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  //print_transaction_list(trans_list,item);
  /* END OF MODIFICATIONS */
  BUF_clear(&LCM_B);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  TRANS_LIST_INIT(&mk_trans_list, LCM_frq, LCM_Ot[item]-LCM_Os[item]);
  /* END OF MODIFICATIONS */
  for(i=0; i<LCM_BM_MAXITEM; i++) perm[i] = LCM_Op[i];
  QUEUE_FE_LOOP_(LCM_jump, i, ii) LCM_Ofrq[ii] = LCM_th;
  LCM_Ofrq[item] = LCM_th;
  /* LAST TWO ARGUMENTS ADDED FOR TRANSACTION KEEPING */
  LCM_mk_freq_trsact(&TT, T, item, LCM_Eend-1, n+(LCM_Ot[item]-LCM_Os[item]), mask, trans_list, &mk_trans_list);
  /* END OF MODIFICATIONS */
  LCM_Ofrq[item] = 0;

  BUF_reset(&LCM_B);
  for(i=0; i<LCM_BM_MAXITEM; i++) LCM_BM_occurrence_delete(i);
  LCMclosed_BM_occurrence_deliver_first(-1, &TT, &mk_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
  LCMclosed_BM_recursive(LCM_BM_MAXITEM, 0xffffffff, 0xffffffff);
  BUF_clear(&LCM_B);

  if(QUEUE_LENGTH(LCM_jump) == 0) goto END0;
  q = ((QUEUE *)(TT.h))->q;
  /* MODIFICATIONS FOR KEEPING TRACK OF TRANSACTIONS */
  if(ii >= 2 && TT.num>5){
	  TRANS_LIST_INIT(&shrink_trans_list, mk_trans_list.siz1, mk_trans_list.siz2);
	  LCM_shrink(&TT, item, 1, &mk_trans_list, &shrink_trans_list);//LAST TWO ARGUMENTS ADDED FOR TRANSACTION KEEPING
	  TRANS_LIST_END(&mk_trans_list);
  }else{
	  shrink_trans_list = mk_trans_list;
  }
  /* END OF MODIFICATIONS */
  LCM_occ_deliver(&TT, item-1);

  do{
    i = QUEUE_ext_tail_(&LCM_jump);
//    printf ("i=%d(%d) %d %d  :%d,%d\n", i, item, LCM_Ot[item]-LCM_Os[item], LCM_Ot[i]-LCM_Os[i], LCM_Ofrq[i], LCM_Ofrq_[i]);
    /* MODIFICATIONS FOR WY ALGORITHM */
    //WY permutations might cause that some items in LCM_jump.q are no longer frequent/testable, so
    // a check must be added
    if(LCM_Ofrq_[i]>=LCM_th) ii = LCMclosed_iter(&TT, i, item, &shrink_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    //else printf("Item %d pruned since LCM_Ofrq_[%d]=%d, and LCM_th=%d\n",i,i,LCM_Ofrq_[i],LCM_th);
    /* END OF MODIFICATIONS */
    LCM_Ot[i] = LCM_Os[i];
    LCM_Ofrq_[i] = 0;
  }while(LCM_jump.t > LCM_jump.s);

  free2(q);
  ARY_end(&TT);
  TRANS_LIST_END(&shrink_trans_list);
  END0:;
  for(i=0; i<LCM_BM_MAXITEM; i++) LCM_Op[i] = perm[i];
  goto END3;
  END2:;
  for(i=LCM_jump.s; i<LCM_jumpt; i++) LCM_Ofrq[LCM_jump.q[i]] = 0;
  LCM_jump.t = LCM_jump.s;
  END3:;
  LCM_jump.s = js;
  while(it<LCM_itemsett) LCM_delitem();
  return (flag);
}

/***************************************************************/
/* main of LCM ver. 3 */
/*************************************************************************/
void LCMclosed(){
  int i;
  BUF_reset(&LCM_B);
  LCMclosed_BM_occurrence_deliver_first(-1, &LCM_Trsact, &root_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
  for (i=LCM_BM_MAXITEM ; i<LCM_Eend ; i++){
    LCM_Ofrq_[i] = LCM_Ofrq[i];
    LCM_Ofrq[i] = 0;
  }
  LCMclosed_BM_recursive(LCM_BM_MAXITEM, 0xffffffff, 0xffffffff);
  BUF_clear(&LCM_B);

  for (i=LCM_BM_MAXITEM ; i<LCM_Eend ; i++){
    //added check with freq
    if(LCM_Ofrq_[i]>=LCM_th) LCMclosed_iter (&LCM_Trsact, i, -1, &root_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    LCM_Ot[i] = LCM_Os[i];
    LCM_Ofrq_[i] = LCM_Ofrq[i] = 0;
  }
  LCM_iters++;
}




// -------------------- NEW METHODS --------------------


void run_significance(int threshold, int num_signif, double deltap, char* trans_file, char* label_file){
  LCM_problem = LCM_CLOSED;
  LCM_print_flag = 1;
  LCM_th = threshold;
  LCM_sign = 1;
  LCM_init (trans_file);
  perm_init(deltap, num_signif, label_file);

  LCMclosed ();

  LCM_end ();
  ARY_end ( &LCM_Trsact );
  transaction_keeping_end();
  perm_end();
}

void run_perm(int threshold, char* trans_file, char* label_file){
  LCM_init (trans_file);
  LCM_th = threshold;
  if(LCM_patterns_with_supp == 0) perm_init(0.0, 0, label_file);
  else perm_reset();

  LCMclosed();

  for(int i=N_BINS-2; i>=0; i--) {
    LCM_r[i] += LCM_r[i+1];
    for(int j=0; j<JP; j++){
      LCM_r_perm[i][j] += LCM_r_perm[i+1][j];
    }
  }
  int* tmp_r = (int*)malloc(JP*sizeof(int));
  for(int i=0; i<N_BINS; i++) {
    memcpy(tmp_r, LCM_r_perm[i], JP*sizeof(int));
    qsort_int(tmp_r, JP);
    LCM_rbeta[i] = tmp_r[(int)(JP*0.95)];
  }
  free(tmp_r);


  LCM_end();
  ARY_end(&LCM_Trsact);
  transaction_keeping_end();
}




int* get_LCM_rbeta(){
  return LCM_rbeta;
}
int* get_LCM_r(){
  return LCM_r;
}
int** get_LCM_r_perm(){
  return LCM_r_perm;
}

int get_LCM_trsact_num(){
  return LCM_trsact_num;
}
int get_LCM_minority_num(char* label_file){
  if(label_file){
    get_N_n(label_file);
  }
  return n;
}
int* get_LCM_patterns_with_supp(){
  return LCM_patterns_with_supp;
}
char* get_LCM_testable_patterns(){
  return LCM_testable_patterns;
}






int get_initial_support(double ub, double target){
	int first = 0, last = N;
	int mid;
	while (first < last) {
		mid = ((last + first) >> 1);
		if (psi[mid]*ub > target)
			first = mid + 1;
		else
			last = mid;
	}
	return first;
}


// --------------- PUBLIC METHODS --------------------
double lambda = 2.0;
int s_start = 1;

void set_lambda(double _lambda){
  lambda = _lambda;
}

void set_seed(int _seed){
  //srand(seed);
  seed = _seed;
}

void set_print_patterns(int flag){
    print_patterns = flag;
}

void set_num_perms(int jp){
    JP = jp;
}

void set_zeta(double zeta){
    lam_bins = 1.0+zeta;
    N_BINS = (int)ceil(log(1e64)/log(lam_bins));
}

void no_sstart(){
  s_start = 0;
}

void run_fastYB(char* trans_file, char* label_file, double alpha){
  tic = second();
  get_N_n(label_file);
  loggamma_init(); psi_init();

  LCM_problem = LCM_CLOSED;
  LCM_print_flag = 0;
  LCM_sign = 0;

  int last_i = N_BINS-1;
  int supp = get_initial_support(1.0, alpha);
  if(!s_start) supp = N/4;
  int besti = last_i;
  double bestfdr = 0;
  int bestsz = 0;
  int besttest = 0;
  free(loggamma); free(psi);
  int cont = 1;
  while(cont && supp>0){
    //printf("Running %d\n", supp); fflush(stdout);
    run_perm(supp, trans_file, label_file);
    //printf("Time %f\n", second()-tic); fflush(stdout);
    int* pws = get_LCM_patterns_with_supp();
    int m_sigma = pws[n];
    for(int i=n-1; i>=0; i--){
      m_sigma += pws[i];
    }

    int* r = get_LCM_r();
    int* rb = get_LCM_rbeta();
    int** rperm = get_LCM_r_perm();

    for(int i=last_i; i>=0; i--){
      //printf("%d %e\n", i, exp(-(i+0.0)*log(lam_bins)));fflush(stdout);
      if( exp(-(i+0.0)*log(lam_bins)) >= psi[supp] ){
        break; // not mined testable patterns at this level
      }
      last_i = i-1;

      double deltap = exp(-(i+0.0)*log(lam_bins));
      double fdr = 0.0, fdr1 = 0.0, fdr2 = 0.0;
      // point estimator
      if(r[i] - rb[i] >=  deltap*m_sigma){ 
        for(int j=0; j<JP; j++){
          fdr1 += ((rperm[i][j]+0.0)/(rperm[i][j]+ r[i] - deltap*m_sigma))/JP;
        }
      }
      else{
        for(int j=0; j<JP; j++){
          fdr1 += (rperm[i][j] > 0 ? 1.0 : 0.0)/JP;
        }
      }
      // upper limit estimator
      if(r[i] - rb[i] > 0){ 
        for(int j=0; j<JP; j++){
          fdr2 += ((rperm[i][j]+0.0)/(rperm[i][j]+ r[i] - rb[i]))/JP;
        }
      }
      else{
        for(int j=0; j<JP; j++){
          fdr2 += (rperm[i][j] > 0 ? 1.0 : 0.0)/JP;
        }
      }
      fdr = fdr1;
      if(fdr2 > fdr1) fdr = fdr2; // get the most conservative estimate

      fdr = fdr/0.95;
      //printf("%d %e\n", i, fdr);fflush(stdout);
      if(fdr <= alpha){
        besti = i;
        bestfdr = fdr;
        bestsz = r[i];
        besttest = m_sigma;
      }
      else{
        cont = 0;
        break;
      }
      //printf("m_sigma %d  ri %d  fdr %f \n", m_sigma, r[i], fdr); fflush(stdout);
    }

    int oldsupp = supp;
    if(cont) supp = floor((supp+0.0)/lambda);
    if(supp == 0 && oldsupp > 1){
      supp = 1;
    }
    if(supp == 0){
      cont = 0;
      break;
    }
  }

  perm_end();

  if(print_patterns) run_significance(supp, 1, exp(-(besti+0.0)*log(lam_bins)), trans_file, label_file);

  //printf("Alpha: %.3f Significant: %d Estimated FDR: %f LastSupp: %d Testable at LastSupp: %d Time %f Peak memory %d\n", alpha, bestsz, bestfdr, supp, besttest, second()-tic, measurePeakMemory());
  printf("%.3f, %d, %f, %d, %d, %f, %d\n", alpha, bestsz, bestfdr, supp, besttest, second()-tic, measurePeakMemory());

}






#endif
