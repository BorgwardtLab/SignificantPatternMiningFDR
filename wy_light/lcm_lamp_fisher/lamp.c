 #ifndef _lamp_c_
#define _lamp_c_

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
// Target FWER
double alpha;
// Number of testable intervals
long long m_testable;
// A (N+1)-dimensional vector such that freq_cnt[j] = #intervals with x_{i}=j processed so far
long long *freq_cnt;

// Region thresholds: Sigma_k = [sl1,sl2] U [N-sl2,N-sl1]
int sl1, sl2;
// Current P-value threshold
double delta;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
int doublecomp(const void*,const void*);
void get_N_n(char *);

// Profiling variables
long long effective_total_dataset_frq;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void lamp_init(double target_fwer, char *labels_file){
	int j; //Loop variable

	get_N_n(labels_file);
	// Ensure class 1 is the minority class
	if(n > (N/2)) n = N-n;

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;
	// And initialise some others
	sl1 = 1; sl2 = N_over_2;
	flag = 1;
	m_testable = 0;
	delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$

	// Allocate memory for freq_cnt, giving an error if it fails
	freq_cnt = (long long *)calloc(N+1, sizeof(long long));
	if(!freq_cnt){
		fprintf(stderr,"Error in function lamp_init: couldn't allocate memory for array freq_cnt\n");
		exit(1);
	}

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	effective_total_dataset_frq = 0; //Init profiling variables
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
void lamp_end(){
	int j, idx_max;
	double delta_corrected;
	delta_corrected = alpha/m_testable;
	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta_corrected);
	fprintf(results_file,"\t Final LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Testable region: [%d,%d] U [%d,%d]\n",sl1,sl2,N-sl2,N-sl1);
	fprintf(results_file,"\t Final P-value lower bound: %e\n",delta);
	fprintf(results_file,"\t Number of testable patterns at final P-value lower bound: %lld\n",m_testable);

	printf("Corrected significance threshold: %e\n",delta_corrected);
	printf("Final LCM support: %d\n",LCM_th);
	printf("Number of testable patterns at final P-value lower bound: %lld\n",m_testable);

	// Free allocated memory
	free(loggamma);
	free(psi);
	free(freq_cnt);

	// Close results file
	fclose(results_file);
}

/* --------------------------------CORE FUNCTIONS------------------------------------ */

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
 *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
 *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
 *    change happened the last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and delta accordingly
 * 3) If sl1 has been modified, then the support of LCM has to be modified
 * 4) Since the temptative corrected significance threshold delta has changed, the FWER needs to be recomputed
 * */
void decrease_threshold(){
	int j; //Loop iterator
	int false_positives; //Number of false positives (a false positive occurs if min_pval[j] <= delta)
	// Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
	if(flag){
		// Update number of testable itemsets
		m_testable -= freq_cnt[sl1]; m_testable -= freq_cnt[N-sl1];
		sl1++; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if(psi[sl1] >= psi[sl2]) delta = psi[sl1];
		else{ delta = psi[sl2]; flag = 0; }
		//Update LCM minimum support
		LCM_th = sl1;
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		if(sl2==(N-sl2)) m_testable -= freq_cnt[sl2];//(beware of case sl2==su1 since it could lead to discounting the same thing twice!)
		else {m_testable -= freq_cnt[sl2]; m_testable -= freq_cnt[N-sl2];}
		sl2--; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; flag = 1; }
		else delta = psi[sl2];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
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

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	effective_total_dataset_frq += x; // Update profiling variable

	// Process testable pattern by increasing counters
	freq_cnt[x]++; m_testable++;

	/* Finally, check if the FWER upper bound constraint is still satisfied, if not decrease threshold */
	while((m_testable*delta) > alpha) {
		//printf("Threshold change BM\n");
		decrease_threshold();
		// Correct possible corruption of LCM data structures due to unexpected change in minimum support
		for(i=0; i<item; i++){
			//printf("Item %d, Frq %d, Current_th %d\n",i,LCM_Ofrq[i],LCM_th);
			if(LCM_Ofrq[i]==(LCM_th-1)){
				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",i,item,LCM_th);
				LCM_BM_occurrence_delete(i);
				*mask &= ~BITMASK_1[i];
				//printf("Problem fixed!\n");
			}
		}
	}
}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	effective_total_dataset_frq += x; // Update profiling variable

	// Process testable pattern by increasing counters
	freq_cnt[x]++; m_testable++;

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while((m_testable*delta) > alpha) {
		//printf("threshold change 0\n");
		decrease_threshold();
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

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	effective_total_dataset_frq += x; // Update profiling variable

	// Sanity-check (this one is more complicated due to the way the transactions are stored)
	aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++) aux++;
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);

	// Process testable pattern by increasing counters
	freq_cnt[x]++; m_testable++;

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while((m_testable*delta) > alpha) {
		//printf("threshold change ary\n");
		decrease_threshold();
		// Correct possible corruption of LCM data structures due to unexpected change in minimum support
		for(j=0; j<LCM_BM_MAXITEM; j++){
			//printf("Item %d, Frq %d, Current_th %d\n",j,LCM_Ofrq[j],LCM_th);
			if(LCM_Ofrq[j]==(LCM_th-1)){
				//printf("Bucket of item %d corrupted after TH change! Parent %d. Current Th%d.\n",j,item,LCM_th);
				LCM_BM_occurrence_delete(j);
				*mask &= ~BITMASK_1[j];
				//printf("Problem fixed!\n");
			}
		}
	}
}

/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
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

#endif
