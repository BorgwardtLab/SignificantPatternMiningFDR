#ifndef _fastyb_c_
#define _fastyb_c_

/* LIBRARY INCLUDES */
#include<math.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"transaction_keeping.c"
#include"lcm_var.c"

/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of the buffer to read chars from file

/* GLOBAL VARIABLES */
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;

// Number of non-empty transaction
int Neff;
// Original vector of labels (dimension # non-empty transactions)
char *labels;
char **labels_perm;

// Corrected significance threshold
double delta;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;

// Output files
//FILE *significant_itemsets_output_file, *pvalues_output_file;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
void pvals_init();
void get_N_n(char *);
void read_labels_file(char *, char*);
extern void LCMFREQ_output_itemset(int *);
static int rand_int(int x);
void randperm(char *buffer, char *src, int l);
double fisher_pval(int a, int x);


double** pre_pval;
int n_pre = -1;
int print_patterns;
/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the code
 * Input arguments are self-explanatory
 * */
void perm_init(double sig_th, int testable_size, char *labels_file){
	// PAOLO
    malloc2(LCM_patterns_with_supp, int, LCM_trsact_num+1, "patterns_with_supp");
    memset(LCM_patterns_with_supp, 0, (LCM_trsact_num+1)*sizeof(int)/sizeof(char));
    LCM_testable_size = testable_size+1;
    if(LCM_print_flag){
        malloc2(LCM_testable_patterns, char, MAX_ITEMSET_STR_LEN, "LCM_testable_patterns");
        memset(LCM_testable_patterns, 0, MAX_ITEMSET_STR_LEN);
    }


	int j; //Loop variable
	char *labels_buffer;

	get_N_n(labels_file);

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	delta = sig_th;

	// Allocate memory for the buffer containing the class labels, giving an error if it fails
	labels_buffer = (char *)malloc(N*sizeof(char));


	/* Allocate memory for the vector of class labels, with labels of empty transactions removed */
	Neff = root_trans_list.siz1;
	labels = (char *)malloc(Neff*sizeof(char));

	// Read file containing class labels and store them in array labels, taking care of removing labels
	// associated with empty transactions
	read_labels_file(labels_file,labels_buffer);
	// Ensure class 1 is the minority class
	if(n > (N/2)){
		for(j=0; j<N; j++) labels_buffer[j] = !labels_buffer[j];
		n = N-n;
	}

	for(j=0;j<Neff;j++) labels[j] = labels_buffer[non_empty_trans_idx[j]];



	srand(seed); // PAOLO change this
	char *perm_buffer = (char *)malloc(N*sizeof(char));;
	labels_perm = (char **)malloc(Neff*sizeof(char *));
	// Now allocate memory for a contiguous block of J*(# non-empty transactions) chars
	labels_perm[0] = (char *)malloc(((long long)JP)*Neff*sizeof(char));
	// And make each row pointer point to the appropriate position (labels_perm[0] is
	// already correctly set)
	for(int i=1;i<Neff;i++) labels_perm[i] = labels_perm[0] + ((long long)i)*JP;

	for(j=0;j<JP;j++) {
		randperm(perm_buffer, labels, N);
		// Dump contents of buffer into destination, skipping values corresponding to empty observations/transactions
		for(int i=0;i<Neff;i++) labels_perm[i][j] = perm_buffer[non_empty_trans_idx[i]];
	}
	free(perm_buffer);

	LCM_r_perm = (int**)malloc(N_BINS*sizeof(int*));
	//LCM_r_perm[0] = (int*)calloc(N_BINS*JP, sizeof(int));
	malloc2(LCM_r_perm[0], int, N_BINS*JP, "LCM_r_perm[0]"); memset(LCM_r_perm[0], 0, N_BINS*JP*sizeof(int));
	for(int i=1; i<N_BINS; i++) LCM_r_perm[i] = LCM_r_perm[0] + ((long long)i)*JP;

	LCM_rbeta = (int*)calloc(N_BINS, sizeof(int));
	LCM_r = (int*)calloc(N_BINS, sizeof(int));

	tmp_ap = (int*)calloc(JP, sizeof(int));



	free(labels_buffer);
	// The array containing the indices of all non-empty transactions is no longer needed
	free(non_empty_trans_idx);

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	//precompute pvalues
	n_pre = MIN(Neff/3, 10000); // limit size to ~1gb
	pvals_init();

}

void perm_reset(){
	memset(LCM_patterns_with_supp, 0, (LCM_trsact_num+1)*sizeof(int)/sizeof(char));
	memset(LCM_r_perm[0], 0, N_BINS*JP*sizeof(int));
	memset(LCM_rbeta, 0, N_BINS*sizeof(int));
	memset(LCM_r, 0, N_BINS*sizeof(int));
}

void perm_end(){
	free(loggamma); free(psi);
	free(pre_pval[0]); free(pre_pval);

	free(tmp_ap);
	free(LCM_r);
	free(LCM_rbeta);
	free(LCM_r_perm[0]);
	free(LCM_r_perm);

	free(labels_perm[0]);
	free(labels_perm);

	free(labels);
    if(LCM_print_flag) free(LCM_testable_patterns);
}




/* Precompute values of log(x!) storing them in the array loggamma */
void loggamma_init(){
	int x;
	// Allocate memory for log-gamma cache, raising error if it fails
	loggamma = (double *)malloc((N+1)*sizeof(double));
	if(!loggamma){
		fprintf(stderr,"Error in function loggamma_init: couldn't allocate memory for array loggamma\n");
		exit(1);
	}
	// Initialise cache with appropriate values
	for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
	// Initialise log_inv_binom_N_n
	log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}

/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	int x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_and_xi1_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	//In this range, the recursion $\psi(x)$=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can be seen to be correct
	for(x=1; x<=n; x++) psi[x] = (((double)(n-(x-1)))/(N-(x-1)))*psi[x-1];

	// Now, start computing xi1 in the range [N-N_over_2,N] using another recursion, this time
	// starting in N
	// Note that we don't need to store all values, since this will be used only to initialise
	// psi[N_over_2]
	x_init = N-N_over_2;
	xi1 = 1;
	//In this range, the recursion $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
	for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n))/(x+1))*xi1;

	// Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get $\psi(N_over_2)$=psi[N_over_2] using the
	// same recursion if N is odd, or simply copy the value of xi1[0] since $\xi_{1}(N-N_over_2)=\psi(N_over_2)$
	// if N is even
	if (N % 2) psi[N_over_2] = (((double)(x_init-n))/x_init)*xi1;
	else psi[N_over_2] = xi1;

	// Now, using psi[N_over_2] compute the right side of "the W", i.e. the range [n+1,N_over_2]
	// using the same recursion as for $\xi_{1}$
	for(x=(N_over_2-1); x > n; x--) psi[x] = (((double)((x+1)-n))/(x+1))*psi[x+1];

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=x_init; x<=N; x++) psi[x] = psi[N-x];

	// Correct minimum attainable P-value in some edge-cases
	if((N % 2)==0){
		if (n == (N/2)) for(x=1; x<N; x++) psi[x] *= 2;
		else psi[N/2] *= 2;
	}
}

void precompute_pvals(double* pre_pval, int x){
	int a_min, a_max, k;
	double p_left, p_right, pval;
	double pre_comp_xterms;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)

	pval = 0; //Accumulate probabilities in this variable
	while(a_min<a_max){
		p_left = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
		p_right = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
		if(p_left == p_right) {
			pval += (p_left+p_right);
			pre_pval[a_min] = pval; pre_pval[a_max] = pval;
			//printf("%d %d  %f\n", x, a_min, log(pval));
			a_min++; a_max--;
		}
		else if(p_left < p_right){
			pval += p_left;
			pre_pval[a_min] = pval;
			a_min++;
		}
		else{
			pval += p_right;
			pre_pval[a_max] = pval;
			a_max--;
		}
	}
	if(a_min == a_max) pre_pval[a_max] = 1.0;
}

void pvals_init(){
	pre_pval = (double**)malloc(n_pre*sizeof(double*));
	pre_pval[0] = (double*)malloc(n_pre*n_pre*sizeof(double));
	for(int i=1; i<n_pre; i++) pre_pval[i] = pre_pval[0] + ((long long)i)*n_pre;

	for(int x=0; x<n_pre; x++){
		int a_min, a_max, k;
		double p_left, p_right, pval;
		double pre_comp_xterms;

		// Compute the contribution of all terms depending on x but not on a
		pre_comp_xterms = loggamma[x] + loggamma[N-x];
		a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
		a_max = (x > n) ? n : x;//min(x,n)

		pval = 0; //Accumulate probabilities in this variable
		while(a_min<a_max){
			p_left = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
			p_right = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
			if(p_left == p_right) {
				pval += (p_left+p_right);
				pre_pval[x][a_min] = pval; pre_pval[x][a_max] = pval;
				//printf("%d %d  %f\n", x, a_min, log(pval));
				a_min++; a_max--;
			}
			else if(p_left < p_right){
				pval += p_left;
				pre_pval[x][a_min] = pval;
				a_min++;
			}
			else{
				pval += p_right;
				pre_pval[x][a_max] = pval;
				a_max--;
			}
		}
		if(a_min == a_max) pre_pval[x][a_max] = 1.0;

	}
}



/* --------------------------------FUNCTIONS TO EVALUATE FISHER'S EXACT TEST P-VALUES ------------------------------------ */

/* Evaluate Fisher's exact test on a table with margins x, n and N and cell count a. Note that n and N are defined as global variables.
 * The p-value is defined as a two-tailed p-value which adds up the probabilities of all tables less or equally likely to occur than the
 * one we observed
 */
double fisher_pval_pre(int a, int x){
	if(x < n_pre) return pre_pval[x][a];
	else return fisher_pval(a, x);
}

double fisher_pval(int a, int x){
	int a_min, a_max, k;
	double p_left, p_right, pval;
	double pre_comp_xterms;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)

	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. As soon as the "accepted" value is located
	// in index a, we know that we have already explored all values of the hypergeometric probability mass whose probabilities are smaller or equal
	// than the probability of a. The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0; //Accumulate probabilities in this variable
	while(a_min<a_max){
		p_left = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
		p_right = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
		if(p_left == p_right) {
			pval += (p_left+p_right);
			if((a==a_min) || (a==a_max)) return pval;
			a_min++; a_max--;
		}
		else if(p_left < p_right){
			pval += p_left;
			if(a==a_min) return pval;
			a_min++;
		}
		else{
			pval += p_right;
			if(a==a_max) return pval;
			a_max--;
		}
	}
	// If we get to this part of the code, it means is the mode of the distribution and, therefore, its associated p-value is 1
	return 1;
}

/* -------------------FUNCTIONS TO PROCESS A NEWLY FOUND TESTABLE HYPOTHESIS-------------------------------------- */

/* This code contains 3 difference functions to process newly found hypotheses. All of them are virtually identical
 * and the only thing which differs is the way the function receives the list of observations (transactions) for
 * which the hypothesis has X=1.
 * LCM has a complex structure, with frequent itemsets being found at 4 different locations in the source code
 * and under 3 different circumstances. Therefore it was needed to introduce differentiation in the way the transaction
 * list is fed to the "solution processing functions" in order to keep the "transaction keeping" overhead minimal.
 *
 * To reuse this code for problems other than frequent itemset mining, the only thing that needs to be modified
 * is the line which computes the cell counts, for example, the following line in bm_process_solution:
 * 		for(i=0; i<current_trans.siz; i++) a += labels_perm[j][current_trans.list[i]];
 * 	There, current_trans.siz is the number of transactions for which the hypothesis has X=1, i.e. the margin x
 * 	of the 2x2 contingency table (note in this case it is redundant with the input argument x of the function)
 * 	Similarly, current_trans.list[i] with i ranging from 0 to (x-1) has the list of indices of the observations
 * 	for which X=1.
 * 	Simply changing those two parameters accordingly allows code reuse.
 * */

/* Process a solution involving the bitmap represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
void bm_process_solution_ap(int x, int item, int *mask){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<current_trans.siz; i++) a += labels[current_trans.list[i]];
	LCM_ap = a;
}
void bm_process_solution_ap_perm(int x, int item, int *mask){
	int i,j;//Loop iterators
	//int* ap = (int*)calloc(JP, sizeof(int));
	memset(tmp_ap, 0, JP*sizeof(int));
	int a = 0;
	char *labels_perm_aux;

	for(i=0; i<current_trans.siz; i++){
		a += labels[current_trans.list[i]];
		labels_perm_aux = labels_perm[current_trans.list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
		for(j=0; !LCM_print_flag & (j<JP); j++){
      tmp_ap[j] += labels_perm_aux[j];
    }
	}
	int cnti=0;
	double* fisher_pvals_pre = (double*)calloc(x+1, sizeof(double));
	if(x < n_pre){
		memcpy(fisher_pvals_pre, pre_pval[x], (x+1)*sizeof(double));
	}
	else{
		precompute_pvals(fisher_pvals_pre, x);
	}
	for(j=0; !LCM_print_flag & (j<JP); j++){
		double pv = fisher_pvals_pre[tmp_ap[j]];//fisher_pval_pre(tmp_ap[j], x);
		int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
		LCM_r_perm[logpv][j]++;
	}
	//printf("bm %d\n", cnti);
	//free(ap);
	LCM_ap = a;

	double pv = fisher_pvals_pre[a];//fisher_pval_pre(a, x);
	int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
	LCM_r[logpv]++;
	free(fisher_pvals_pre);
}


/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0_ap(int x){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<bm_trans_list[1].siz; i++) a += labels[bm_trans_list[1].list[i]];
	LCM_ap = a;
}
void process_solution0_ap_perm(int x){
	int i,j;//Loop iterators
	//int* ap = (int*)calloc(JP, sizeof(int));
	memset(tmp_ap, 0, JP*sizeof(int));
	int a = 0;
	char *labels_perm_aux;

	for(i=0; i<bm_trans_list[1].siz; i++){
		a += labels[bm_trans_list[1].list[i]];
		labels_perm_aux = labels_perm[bm_trans_list[1].list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
		for(j=0; !LCM_print_flag & (j<JP); j++){
      tmp_ap[j] += labels_perm_aux[j];
    }
	}
	int cnti=0;
	double* fisher_pvals_pre = (double*)calloc(x+1, sizeof(double));
	if(x < n_pre){
		memcpy(fisher_pvals_pre, pre_pval[x], (x+1)*sizeof(double));
	}
	else{
		precompute_pvals(fisher_pvals_pre, x);
	}
	for(j=0; !LCM_print_flag & (j<JP); j++){
		double pv = fisher_pvals_pre[tmp_ap[j]];//fisher_pval_pre(tmp_ap[j], x);
		int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
		LCM_r_perm[logpv][j]++;
	}
	//printf("s0 %d\n", cnti);
	//free(ap);
	LCM_ap = a;
	double pv = fisher_pvals_pre[a];//fisher_pval_pre(a, x);
	int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
	LCM_r[logpv]++;
	free(fisher_pvals_pre);
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution_ap(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	int a; //Cell count of current itemset
	// Compute the cell-count corresponding to the current itemset, plus sanity-check
	a = 0; aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
			a += labels[*ptr];
			aux++;
		}
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);
	LCM_ap = a;
}
void ary_process_solution_ap_perm(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	int a; //Cell count of current itemset
	//int* ap = (int*)calloc(JP, sizeof(int));
	memset(tmp_ap, 0, JP*sizeof(int));
	char *labels_perm_aux;

	// Compute the cell-count corresponding to the current itemset, plus sanity-check
	a = 0; aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
			labels_perm_aux = labels_perm[*ptr];
			for(j=0; !LCM_print_flag & (j<JP); j++){
        tmp_ap[j] += labels_perm_aux[j];
      }

			a += labels[*ptr];
			aux++;
		}
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);
	int cnti=0;
	double* fisher_pvals_pre = (double*)calloc(x+1, sizeof(double));
	if(x < n_pre){
		memcpy(fisher_pvals_pre, pre_pval[x], (x+1)*sizeof(double));
	}
	else{
		precompute_pvals(fisher_pvals_pre, x);
	}
	for(j=0; !LCM_print_flag & (j<JP); j++){
		double pv = fisher_pvals_pre[tmp_ap[j]];//fisher_pval_pre(tmp_ap[j], x);
		//printf("%f ", log(pv));
		int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
		LCM_r_perm[logpv][j]++;
	}
	//printf("bm %d\n", cnti);
	//free(ap);
	LCM_ap = a;
	double pv = fisher_pvals_pre[a];//fisher_pval_pre(a, x);
	int logpv = floor(-log(pv + 1e-60)/log(lam_bins) +1e-10);
	LCM_r[logpv]++;
	free(fisher_pvals_pre);
}





/* Do a first scan of the file containing the class labels to compute the total number of observations, N,
 * and the total number of observations in the positive class, n
 * */
void get_N_n(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	N = 0; n = 0;

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function get_N_n when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			N++;
			if(char_to_int[*read_buf_aux]) n++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

void read_labels_file(char *labels_file, char *labels_buffer){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	char *labels_aux = labels_buffer;//Auxiliary pointer to array labels for increments

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function read_labels_file when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			*labels_aux++ = char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}




static int rand_int(int x){
	int rnd;
	int limit = RAND_MAX - RAND_MAX % x;

	do{
		rnd = rand();
	}while(rnd >= limit);
	return rnd % x;
}

void randperm(char *buffer, char *src, int l){
	int i,j; // Variables for looping and swapping
	char tmp; // Temp int for swapping

	// First of all, copy the original array in the buffer
	for(i=0;i<l;i++) buffer[i] = src[i];

	// Fisher-Yates algorithm
	for(i = l-1; i > 0; i--){
		// Sample a random integer in [0,i]
		j = rand_int(i + 1);
		// Swap dest[j] and dest[i]
		tmp = buffer[j];
		buffer[j] = buffer[i];
		buffer[i] = tmp;

	}
}

#endif
