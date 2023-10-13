#ifndef _transaction_keeping_c_
#define _transaction_keeping_c_

/* INCLUDE DEPENDENCIES ON LCM CODE */
//Need access to some internal LCM variables
#include"transaction_keeping.h"
#include"var_declare.h"

/* MACROS */
// Compute rounded up integer of quotient a/b
#define CEIL(a, b) (((a) / (b)) + (((a) % (b)) > 0 ? 1 : 0))

/* GLOBAL VARIABLES */
TRANS_LIST root_trans_list;
BM_TRANS_LIST *bm_trans_list;
BM_TRANS_LIST current_trans;
int **shrink_workspace1;
int *workspace1;
int *bitmap_item_frq;
int bm_trans_list_nodes;
int *non_empty_trans_idx;
//int print_counter = 0;

/* FUNCTION DECLARATIONS */
void remove_empty_transactions(char *);

/* AUXILIARY FUNCTIONS */

void TRANS_LIST_INIT(TRANS_LIST *L,int siz1,int siz2){
	L->siz1 = siz1;
	L->siz2 = siz2;
	L->list = (int *)malloc(siz1*sizeof(int));
	if(!L->list){
		fprintf(stderr,"Error in function TRANS_LIST_INIT: couldn't allocate memory for array L->list\n");
		exit(1);
	}
	L->ptr = (int **)malloc(siz2*sizeof(int *));
	if(!L->ptr){
		fprintf(stderr,"Error in function TRANS_LIST_INIT: couldn't allocate memory for array L->ptr\n");
		exit(1);
	}
}

void TRANS_LIST_END(TRANS_LIST *L){
	free(L->list);
	free(L->ptr);
}

void BM_TRANS_LIST_INIT(int n_items){
	int i,n_nodes,mem_siz;
	int l1,l2;
	//Compute 2^LCM_BM_MAXITEM
	for(i=0,n_nodes=1;i<n_items;i++,n_nodes*=2);
	//Allocate memory for BM_TRANS_LIST structure
	bm_trans_list = (BM_TRANS_LIST *)malloc(n_nodes*sizeof(BM_TRANS_LIST));
	if(!bm_trans_list){
		fprintf(stderr,"Error in function BM_TRANS_LIST_INIT: couldn't allocate memory for array bm_trans_list\n");
		exit(1);
	}
	//Initialize lists inside each node of BM_TRANS_LIST
	//The (i+1)-th most frequent item has 2^i nodes, hence each node gets a share
	//of frq[i]/2^i *sizeof(int) bytes
	for(i=0,l1=1,l2=2;i<n_items;i++,l2*=2){
		mem_siz = CEIL(bitmap_item_frq[i],l1);
		//mem_siz = bitmap_item_frq[i];
		for(;l1<l2;l1++){
			bm_trans_list[l1].list = (int *)malloc(mem_siz*sizeof(int));
			if(!bm_trans_list[l1].list){
				fprintf(stderr,"Error in function BM_TRANS_LIST_INIT: couldn't allocate memory for array bm_trans_list[l1].list\n");
				exit(1);
			}
			bm_trans_list[l1].max_siz = mem_siz;
			bm_trans_list[l1].siz = 0;
		}
	}
	bm_trans_list_nodes = n_nodes;
	// Initialize structure to keep current transaction
	current_trans.list = (int *)malloc(bitmap_item_frq[0]*sizeof(int));
	if(!current_trans.list){
		fprintf(stderr,"Error in function BM_TRANS_LIST_INIT: couldn't allocate memory for array current_trans.list\n");
		exit(1);
	}
	current_trans.max_siz = bitmap_item_frq[0];
	current_trans.siz = 0;
}

void BM_TRANS_LIST_END(){
	int i;
	// Note that bm_trans_list[0] is in fact never initialized
	for(i=1;i<bm_trans_list_nodes;i++) free(bm_trans_list[i].list);
	free(bm_trans_list);
}

void BM_TRANS_LIST_INSERT(int p, int *src, int siz){
	int new_size;
	// If p==0, the transaction does not belong in the CPT
	if(p==0) return;
	// Check if current list size if big enough to fit new data
	// If not, allocate double the current memory size
	new_size = bm_trans_list[p].siz + siz;
	if(new_size > bm_trans_list[p].max_siz){
		new_size = (new_size > 2*bm_trans_list[p].max_siz) ? new_size : 2*bm_trans_list[p].max_siz;
		bm_trans_list[p].list = (int *)realloc(bm_trans_list[p].list,new_size*sizeof(int));
		if(!bm_trans_list[p].list){
			fprintf(stderr,"Error in function BM_TRANS_LIST_INSERT: couldn't reallocate memory for array bm_trans_list[p].list\n");
			exit(1);
		}
		bm_trans_list[p].max_siz = new_size;
	}
	memcpy(bm_trans_list[p].list + bm_trans_list[p].siz,src,siz*sizeof(int));
	bm_trans_list[p].siz += siz;
}

void BM_TRANS_LIST_EMPTY(int p){
	bm_trans_list[p].siz = 0;
}

void BM_CURRENT_TRANS_INSERT(int p){
	// This list can never overflow its initially allocated size, equal to the frequency of the
	// most frequent item
	memcpy(current_trans.list + current_trans.siz,bm_trans_list[p].list,bm_trans_list[p].siz*sizeof(int));
	current_trans.siz += bm_trans_list[p].siz;
}

void BM_CURRENT_TRANS_EMPTY(){
	current_trans.siz = 0;
}

/* INITIALIZATION AND TERMINATION FUNCTIONS */

void transaction_keeping_init(char *trans_file){
	int i;
	non_empty_trans_idx = (int *)malloc(LCM_Trsact.num*sizeof(int));
	if(!non_empty_trans_idx){
		fprintf(stderr,"Error in function transaction_keeping_init: couldn't allocate memory for array non_empty_trans_idx\n");
		exit(1);
	}
	remove_empty_transactions(trans_file);
	TRANS_LIST_INIT(&root_trans_list,LCM_Trsact.num,LCM_Trsact.num);
	for(i=0;i<LCM_Trsact.num;i++){
		root_trans_list.list[i] = i;
		root_trans_list.ptr[i] = &root_trans_list.list[i];
	}
	// Cipher-permutation version
	/*
	for(i=0;i<LCM_Trsact.num;i++){
		root_trans_list.list[i] = non_empty_trans_idx[i];
		root_trans_list.ptr[i] = &root_trans_list.list[i];
	}
	// non_empty_trans_idx no longer necessary
	free(non_empty_trans_idx);
	*/
	workspace1 = (int *)malloc(LCM_Trsact.num*sizeof(int));
	if(!workspace1){
		fprintf(stderr,"Error in function transaction_keeping_init: couldn't allocate memory for array workspace1\n");
		exit(1);
	}
	shrink_workspace1 = (int **)malloc(LCM_Trsact.num*sizeof(int *));
	if(!shrink_workspace1){
		fprintf(stderr,"Error in function transaction_keeping_init: couldn't allocate memory for array shrink_workspace1\n");
		exit(1);
	}
	for(i=0;i<LCM_Trsact.num;i++) shrink_workspace1[i] = ((int *)0);
}

void transaction_keeping_end(){
	TRANS_LIST_END(&root_trans_list);
	BM_TRANS_LIST_END();
	free(bitmap_item_frq);
	free(shrink_workspace1);
	free(workspace1);
	free(current_trans.list);
	//printf("Number of lines outputted: %d\n",print_counter);
}

/* FUNCTION FOR DEALING WITH EMPTY TRANSACTIONS */
void remove_empty_transactions(char *trans_file){
	char c;
	FILE *f_trans = ((FILE*)0);
	int n_lines = 0, n_lines_non_empty = 0;
	int non_empty_line_flag = 0;

	//Try to open file, giving an error message if it fails
	if(!(f_trans = fopen(trans_file,"r"))){
		fprintf(stderr, "Error in function remove_empty_transactions when opening file %s\n",trans_file);
		exit(1);
	}

	// Scan transaction database to find non-empty transactions and store their indices in
	// non_empty_trans_idx
	while((c = fgetc(f_trans)) != EOF){
		if(c == '\n'){
			if(non_empty_line_flag) non_empty_trans_idx[n_lines_non_empty++] = n_lines;
			n_lines++;
			non_empty_line_flag = 0;
		}else{
			non_empty_line_flag = 1;
		}
	}

	// Close file
	fclose(f_trans);
}


/* FUNCTION FOR DISPLAYING TRANSACTIONS OF A NON-BITMAP REPRESENTED ITEMSET */
/*
void print_transaction_list(TRANS_LIST *L,int item){
	int *t, *t_end;
	int *ptr, *end_ptr;
	print_counter++;
	return;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
			if((ptr==(end_ptr-1)) && (t==(t_end-1))) printf("%d",*ptr);
			else printf("%d,",*ptr);
		}
	}
	printf("\n");
}

void print_current_trans(){
	int i;
	print_counter++;
	return;
	for(i=0;i<(current_trans.siz-1);i++) printf("%d,",current_trans.list[i]);
	printf("%d\n",current_trans.list[i]);
}

void print_trans0(){
	int i;
	print_counter++;
	return;
	for(i=0;i<(bm_trans_list[1].siz-1);i++) printf("%d,",bm_trans_list[1].list[i]);
	printf("%d\n",bm_trans_list[1].list[i]);
}
*/
#endif
