#ifndef _enumerate_significant_itemsets_c_
#define _enumerate_significant_itemsets_c_

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
FILE* results_file;
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;

// Number of non-empty transaction
int Neff;
// Original vector of labels (dimension # non-empty transactions)
char *labels;

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
FILE *significant_itemsets_output_file, *pvalues_output_file;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
void get_N_n(char *);
void read_labels_file(char *, char*);
extern void LCMFREQ_output_itemset(int *);
// Profiling variables
long long n_significant_patterns;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the code
 * Input arguments are self-explanatory
 * */
void enum_sig_itemsets_init(double sig_th, char *labels_file){
	int j; //Loop variable
	char *labels_buffer;

	get_N_n(labels_file);

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	delta = sig_th;

	// Allocate memory for the buffer containing the class labels, giving an error if it fails
	labels_buffer = (char *)malloc(N*sizeof(char));
	if(!labels_buffer){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels_buffer\n");
		exit(1);
	}

	/* Allocate memory for the vector of class labels, with labels of empty transactions removed */
	Neff = root_trans_list.siz1;
	labels = (char *)malloc(Neff*sizeof(char));
	if(!labels){
		fprintf(stderr,"Error in function enum_sig_itemsets_init: couldn't allocate memory for array labels\n");
		exit(1);
	}

	// Read file containing class labels and store them in array labels, taking care of removing labels
	// associated with empty transactions
	read_labels_file(labels_file,labels_buffer);
	// Ensure class 1 is the minority class
	if(n > (N/2)){
		for(j=0; j<N; j++) labels_buffer[j] = !labels_buffer[j];
		n = N-n;
	}

	for(j=0;j<Neff;j++) labels[j] = labels_buffer[non_empty_trans_idx[j]];
	free(labels_buffer);
	// The array containing the indices of all non-empty transactions is no longer needed
	free(non_empty_trans_idx);

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	// Initialise profiling variables
	n_significant_patterns = 0;
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

/* Free all allocated memory and give some output for debugging purposes */
void enum_sig_itemsets_end(){
	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta);
	fprintf(results_file,"\t LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Number of significant patterns found: %lld\n",n_significant_patterns);

	printf("Number of significant patterns found: %lld\n",n_significant_patterns);

	// Free allocated memory
	free(loggamma); free(psi);
	free(labels);

	// Close output files
	fclose(results_file);
	fclose(significant_itemsets_output_file);
	fclose(pvalues_output_file);
}

/* --------------------------------FUNCTIONS TO EVALUATE FISHER'S EXACT TEST P-VALUES ------------------------------------ */

/* Evaluate Fisher's exact test on a table with margins x, n and N and cell count a. Note that n and N are defined as global variables.
 * The p-value is defined as a two-tailed p-value which adds up the probabilities of all tables less or equally likely to occur than the
 * one we observed
 */
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
void bm_process_solution(int x, int item, int *mask){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	double pval;//P-value of current itemset

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);
	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;
	
	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<current_trans.siz; i++) a += labels[current_trans.list[i]];
	// Compute the corresponding p-value
	pval = fisher_pval(a,x);
	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	//fprintf_current_itemset(); //remove this
	//printf("%d %d %d %e %f\n", a, x, current_trans.siz, pval, delta);
	if(pval <= delta){
		n_significant_patterns++;
		//fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
		fprintf_current_itemset();
	}


}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators
	int a; //Cell count of current itemset
	double pval;//P-value of current itemset

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);
	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;
	// Compute the cell-count corresponding to the current itemset
	a = 0;
	for(i=0; i<bm_trans_list[1].siz; i++) a += labels[bm_trans_list[1].list[i]];
	// Compute the corresponding p-value
	pval = fisher_pval(a,x);
	//LCMFREQ_output_itemset(LCM_add.q+LCM_add.t); //remove this
	//printf("%d %d %d %e %f\n", a, x, bm_trans_list[1].siz, pval, delta);
	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	if(pval <= delta){
		n_significant_patterns++;
		//fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
		LCMFREQ_output_itemset(LCM_add.q+LCM_add.t);
	}
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	int a; //Cell count of current itemset
	double pval;//P-value of current itemset

	/* First, process the new hypothesis */
	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;
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
	// Compute the corresponding p-value
	pval = fisher_pval(a,x);
	//fprintf_current_itemset(); //remove this
	//printf("%d %d %e %f  %d %d\n", a, x, pval, delta, LCM_Os[item], LCM_Ot[item]);
	// If p-value is significant, write current itemset and the corresponding p-value to the output files
	if(pval <= delta){
		n_significant_patterns++;
		//fprintf(pvalues_output_file,"%d,%d,%.18e\n",a,x,pval);
		fprintf_current_itemset();
	}

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

#endif
