/* Linear time Closed itemset Miner for Frequent Itemset Mining problems */
/* 2004/4/10 Takeaki Uno */
/* This program is available for only academic use.
   Neither commercial use, modification, nor re-distribution is allowed */

#ifndef _lcm_c_
#define _lcm_c_

#include<time.h>
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

/* MODIFICATIONS FOR LAMP ALGORITHIM */
#include"lamp.c"
/* END OF MODIFICATIONS */

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
  for(i=1; i<item; i++)
    if(LCM_Ofrq[i] >= LCM_th) LCMclosed_BM_iter(i, mask, pmask);
    else{
    	if(LCM_Ofrq[i]>0) {
    		//printf("Item %d is infrequent (%d,%d) and will not have LCMclosed_BM_iter,\n",i,LCM_Ofrq[i],LCM_th);
    		LCM_BM_occurrence_delete(i);
    		//printf("Problem fixed!\n");
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

  LCM_solution();

  /* MODIFICATION FOR FAST WY ALGORITHIM */
  if(LCM_frq != current_trans.siz){
	  printf("LCM_frq=%d, current_trans.siz=%d\n",LCM_frq,current_trans.siz);
  }
  bm_process_solution(LCM_frq,item,&mask);
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
    LCM_solution();
    /* MODIFICATIONS FOR WY ALGORITHM */
    ary_process_solution(LCM_frq, trans_list, item, &mask);
    /* END OF MODIFICATIONS */
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
  LCM_solution();
  /* MODIFICATIONS FOR WY ALGORITHM */
  ary_process_solution(LCM_frq, trans_list, item, &mask);
  /* END OF MODIFICATIONS */
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
    LCMclosed_iter (&LCM_Trsact, i, -1, &root_trans_list);//LAST ARGUMENT ADDED FOR TRANSACTION KEEPING
    LCM_Ot[i] = LCM_Os[i];
    LCM_Ofrq_[i] = LCM_Ofrq[i] = 0;
  }
  LCM_iters++;
}

/*************************************************************************/
/*************************************************************************/
int main(int argc, char *argv[]){
  int i;

  /* MODIFICATIONS FOR FAST WY ALGORITHIM */

  // Main input arguments which are not part of LCM
  double sig_th;
  char *class_labels_file;
  char *tmp_filename;

  // Initial time
  t_init = measureTime();

  // Check if input contains all needed arguments
  if (argc != 5){
	  printf("LCM_LAMP_FISHER: output_basefilename target_fwer input_class_labels_file input_transactions_file\n");
	  exit(1);
  }

  // Create output files for results and profiling
  tmp_filename = (char *)malloc((strlen(argv[1])+512)*sizeof(char));
  if(!tmp_filename){
	fprintf(stderr,"Error in function main: couldn't allocate memory for array tmp_filename\n");
	exit(1);
  }
  // Create a file to report runtime information
  strcpy(tmp_filename,argv[1]); strcat(tmp_filename,"_timing.txt");
  if(!(timing_file = fopen(tmp_filename,"w"))){
	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
	exit(1);
  }
  // Create a file to report results
  strcpy(tmp_filename,argv[1]); strcat(tmp_filename,"_results.txt");
  if(!(results_file = fopen(tmp_filename,"w"))){
  	fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
  	exit(1);
  }

  // Free filename holder
  free(tmp_filename);

  // Additional arguments
  sig_th = atof(argv[2]);
  class_labels_file = argv[3];

  // Remove arguments introduced by WY functionality to avoid interference with the rest of LCM's code
  // Only remaining arguments at this point should be transactions_file = argv[1]
  argv[1] = argv[4]; argc = 2;

  // Initialise the support of LCM to 1
  LCM_th = 1;
  /* END OF MODIFICATIONS */
  tic = measureTime();
  LCM_problem = LCM_CLOSED;
  LCM_init(argc, argv);
  toc = measureTime();
  time_LCM_init = toc-tic;

  /* MODIFICATIONS FOR FAST WY ALGORITHIM */

  // Initialize Westfall-Young permutation code
  tic = measureTime();
  lamp_init(sig_th,class_labels_file);
  toc = measureTime();
  time_initialisation_lamp = toc-tic;
  /* END OF MODIFICATIONS */
  tic = measureTime();
  LCMclosed();
  toc = measureTime();
  time_threshold_correction = toc-tic;

  // Main part of the code
  LCM_output();
  LCM_end();
  ARY_end(&LCM_Trsact);

  /* MODIFICATION TO KEEP TRACK OF TRANSACTIONS */
  tic = measureTime();
  transaction_keeping_end();
  /* END OF MODIFICATIONS */
  /* MODIFICATIONS FOR FAST WY ALGORITHIM */
  lamp_end();
  toc = measureTime();
  time_termination = toc-tic;
  // Final time
  t_end = measureTime();
  /* END OF MODIFICATIONS */

  /* MODIFICATIONS FOR CODE PROFILING */
  profileCode();
  /* END OF MODIFICATIONS */

  exit(0);
}


#endif


