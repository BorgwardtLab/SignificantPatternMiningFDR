#ifndef _transaction_keeping_h_
#define _transaction_keeping_h_

/* MAIN STRUCT FOR TRANSACTION KEEPING */

typedef struct{
	int siz1;//Number of original transactions in the list
	int siz2;//Number of merged transactions in the list
	int *list;//List of all original transactions
	int **ptr;//List of pointers such that ptr[i] points to the point of list such that the transactions belonging to the i-th merged transaction begins
}TRANS_LIST;

typedef struct{
	int siz;//Current size
	int max_siz;//Maximum allocated memory
	int *list;//Pointer to contents
}BM_TRANS_LIST;

#endif
