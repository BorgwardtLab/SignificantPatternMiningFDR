#ifndef _wy_c_
#define _wy_c_

/* LIBRARY INCLUDES */
#include<math.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include"var_declare.h"
#include"permutation.c"
#include"transaction_keeping.c"
#include"lcm_var.c"

/* CONSTANT DEFINES */


/* GLOBAL VARIABLES */
FILE* results_file, *minpvals_file;
// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
int N, N_over_2;
// Number of observations in positive class
int n;
// Target FWER
double alpha;
// Current FWER
double FWER;

// Minimum P-value for each permutation
double *min_pval;
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
// Array for storing values of the PDF of the hypergeometric distribution for fast p-value computation
// and eventually the p-values themselves (the array is reused)
double *hypergeom_pvals;

// Cell-count counter
int *a_cnt;

/* FUNCTION DECLARATIONS */
void precompute_pvals(int);
void loggamma_init();
void psi_init();
int doublecomp(const void*,const void*);

// Profiling variables
long long n_pvalues_computed;
long long n_cellcounts_computed;
long long effective_total_dataset_frq;

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void wy_init(double target_fwer){
	int j; //Loop variable

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;
	// And initialise some others
	sl1 = 1; sl2 = N_over_2;
	flag = 1;
	FWER = 0;
	delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$

	// Initialise cache for log(x!) and psi(x)
	loggamma_init();
	psi_init();

	// Allocate memory for minimum p-values, raising error if it fails
	min_pval = (double *)malloc(J*sizeof(double));
	if(!min_pval){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array min_pval\n");
		exit(1);
	}
	// Initialise all p-values to 1
	for(j=0; j<J; j++) min_pval[j] = 1;

	// Allocate memory for the precomputed p-values of the hypergeometric distribution
	// (worst case memory requirement n+1), raising an error if it fails
	hypergeom_pvals = (double *)malloc((n+1)*sizeof(double));
	if(!hypergeom_pvals){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array hypergeom_pvals\n");
		exit(1);
	}

	// Allocate memory for cell counts, raising an error if it fails
	a_cnt = (int *)malloc(J*sizeof(int));
	if(!a_cnt){
		fprintf(stderr,"Error in function wy_init: couldn't allocate memory for array a_cnt\n");
		exit(1);
	}
	for(j=0; j<J; j++) a_cnt[j] = 0;

	n_pvalues_computed = 0; n_cellcounts_computed = 0; effective_total_dataset_frq = 0; //Init profiling variables
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
void wy_end(){
	int j, idx_max;
	double delta_corrected;

	// Sort p-values
	qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_corrected = min_pval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(delta_corrected==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_corrected);
		delta_corrected = min_pval[idx_max];
	}

	// Print results
	fprintf(results_file,"RESULTS\n");
	fprintf(results_file,"\t Corrected significance threshold: %e\n",delta_corrected);
	fprintf(results_file,"\t FWER at corrected significance threshold: %e\n",floor(idx_max+1)/J);
	fprintf(results_file,"\t Final LCM support: %d\n",LCM_th);
	fprintf(results_file,"\t Final P-value lower bound: %e\n",delta);
	fprintf(results_file,"\t FWER at final P-value lower bound: %e\n",FWER);

	fprintf(minpvals_file,"MINIMUM P-VALS (%d PERMUTATIONS)\n",J);
	for(j=0;j<(J-1);j++) fprintf(minpvals_file,"%e,",min_pval[j]);
	fprintf(minpvals_file,"%e\n",min_pval[J-1]);

	// Report the seed for reproducibility
	fprintf(minpvals_file,"\nRANDOM SEED USED FOR KEY GENERATION\n");
	fprintf(minpvals_file,"\t Seed = %u\n",(unsigned)seed);

	// Free allocated memory
	free(loggamma);
	free(psi);
	free(min_pval);
	free(hypergeom_pvals);
	free(a_cnt);

	printf("Corrected significance threshold: %e\n",delta_corrected);
	printf("Final LCM support: %d\n",LCM_th);
	//printf("Number of testable patterns at final P-value lower bound: %lld\n",m_testable);

	// Close files
	fclose(results_file); fclose(minpvals_file);
}

/* -------------------------------FUNCTIONS TO COMPUTE FISHER EXACT TEST P-VALUES----------------------------------- */

/* This function precomputes all Fisher exact test P-values for a contingency table with margins x,n,N that is,
 * all p-values p(a,x,n,N) for a in the range [max(0,n+x-N),min(x,n)]. The results will be stored in the array
 * hypergeom_pvals such that p(a,x,n,N)=hypergeom_pvals[a]. Note that values hypergeom_pvals[a] for a outside
 * [max(0,n+x-N),min(x,n)] are undefined and could contain garbage of previous hypotheses.
 * */
void precompute_pvals(int x){
	double pre_comp_xterms, pval, p_left, p_right;
	int a, a_min, a_max;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)

	// Precompute the hypergeometric PDF in the range of interest
	for(a=a_min; a<=a_max; a++) hypergeom_pvals[a] = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));

	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. When a value is "accepted", we know its
	// respective p-value because due to the way we explore the hypergeometric pdf, there can be no other values larger than it. Therefore, everytime
	// a value is "accepted", we store the pvalue in hypergeom_pvals.
	// The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0;
	while(a_min<a_max){
		p_left = hypergeom_pvals[a_min]; p_right = hypergeom_pvals[a_max];
		if(p_left == p_right) { pval += (p_left+p_right); hypergeom_pvals[a_min++] = pval; hypergeom_pvals[a_max--] = pval; }
		else if(p_left < p_right){ pval += p_left; hypergeom_pvals[a_min++] = pval;}
		else{ pval += p_right; hypergeom_pvals[a_max--] = pval;}
	}
	// In this case a_min=a_max is the mode of the distribution and its p-value is 1 by definition
	if(a_min==a_max) hypergeom_pvals[a_max] = 1;
}


/* --------------------------------CORE FAST WESTFALL-YOUNG PERMUTATION FUNCTIONS------------------------------------ */

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
		sl1++; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]) delta = psi[sl1];
		else{ delta = psi[sl2]; flag = 0; }
		//Update LCM minimum support
		LCM_th = sl1;
		//printf("\n\n\nTHRESHOLD CHANGE!!! NEW THRESHOLD=%d\n\n\n",LCM_th);
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl2--; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; flag = 1; }
		else delta = psi[sl2];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	// Recompute FWER from scratch
	false_positives = 0;
	for(j=0; j<J; j++) false_positives += (min_pval[j]<=delta) ? 1 : 0;
	FWER = ((double)false_positives)/J;
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
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	// Sanity-check
	if (x != current_trans.siz) printf("Error: x = %d, current_trans.siz=%d\n",x,current_trans.siz);

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute cell-counts for all J-permutations
	for(i=0; i<current_trans.siz; i++){
		labels_perm_aux = labels_perm[current_trans.list[i]];//Avoid recomputing labels_perm[current_trans.list[i]] all the time
		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in bm_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
			if( (pval<=delta) && (min_pval[j]>delta)) FWER += ((double)1)/J;
			min_pval[j] = pval;
		}
	}

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
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
	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif
}

/* Process a solution involving the most frequent item (item 0) */
// x = frequency (i.e. number of occurrences) of newly found solution
void process_solution0(int x){
	int i,j;//Loop iterators
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	// Sanity-check
	if (x != bm_trans_list[1].siz) printf("Error: x = %d, bm_trans_list[1].siz=%d\n",x,bm_trans_list[1].siz);

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute cell-counts for all permutations
	for(i=0; i<bm_trans_list[1].siz; i++){
		labels_perm_aux = labels_perm[bm_trans_list[1].list[i]]; //Avoid recomputing labels_perm[bm_trans_list[1].list[i]] all the time
		for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in process_solution0!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
			if( (pval<=delta) && (min_pval[j]>delta)) FWER += ((double)1)/J;
			min_pval[j] = pval;
		}
	}
	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
		//printf("threshold change 0\n");
		decrease_threshold();
	}

	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif
}

/* Process a solution involving the array-list represented itemsets */
// x = frequency (i.e. number of occurrences) of newly found solution
// L = pointer to TRANS_LIST struct keeping track of merged transactions
// item = current node of the tree
void ary_process_solution(int x, TRANS_LIST *L, int item, int *mask){
	int j;//Loop iterator
	int aux; //Auxiliary counter
	int *t, *t_end, *ptr, *end_ptr; //Pointers for iterating on transaction list
	double pval; //Variable to hold p-values
	char *labels_perm_aux; //Auxiliary pointer

	/* First, process the new hypothesis */

	// Minimum attainable P-value for the hypothesis
	double psi_x = psi[x];
	// Check if the newly found solution is in the current testable region Sigma_k
	if(psi_x > delta) return;

	#ifdef PROFILE_MINING
	ticp = measureTime();
	#endif

	// Sanity-check (this one is more complicated due to the way the transactions are stored)
	aux = 0;
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++) aux++;
	}
	if (x != aux) printf("Error: x = %d, trans_size=%d\n",x,aux);

	// Precompute CDF and p-values of hypergeometric distribution with frequency x
	precompute_pvals(x);
	n_pvalues_computed++; //Update profiling variable

	// Compute cell-counts for all permutations
	for(t=LCM_Os[item],t_end=LCM_Ot[item];t<t_end;t++){
		end_ptr = (*t == (L->siz2-1)) ? L->list + L->siz1 : L->ptr[*t+1];
		for(ptr = L->ptr[*t];ptr < end_ptr;ptr++){
			labels_perm_aux = labels_perm[*ptr];
			for(j=0; j<J; j++) a_cnt[j] += labels_perm_aux[j];
		}
	}
	n_cellcounts_computed += J; //Update profiling variable
	effective_total_dataset_frq += x; // Update profiling variable

	// If not, compute permuted P-values for each permutation
	for(j=0; j<J; j++){
		// Fetch the precomputed p-value
		pval = hypergeom_pvals[a_cnt[j]];
		// Sanity-check
		if(pval < 0) printf("Negative P-value detected in ary_process_solution!: j=%d, x=%d, a=%d, pval=%e.\n",j,x,a_cnt[j],pval);
		a_cnt[j] = 0;
		// Check if the newly computed p-value is smaller than current minimum p-value for the
		// permutation
		if(pval < min_pval[j]){
			// Check if the decrease in the current minimum p-value for the j-th permutation
			// causes an increase in the FWER
			if( (pval<=delta) && (min_pval[j]>delta)) FWER += ((double)1)/J;
			min_pval[j] = pval;
		}
	}

	/* Finally, check if the FWER constraint is still satisfied, if not decrease threshold */
	while(FWER > alpha) {
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

	#ifdef PROFILE_MINING
	tocp = measureTime(); time_minpval += tocp-ticp;
	#endif

}

/* AUXILIARY FUNCTIONS */
// Comparison function used by quicksort implementation in C
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

#endif
