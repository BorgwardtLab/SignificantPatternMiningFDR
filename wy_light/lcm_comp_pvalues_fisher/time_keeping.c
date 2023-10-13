#ifndef _time_keeping_c_
#define _time_keeping_c_

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
#include <sys/time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

/* CODE DEPENDENCIES */
#include"var_declare.h"

/* GLOBAL VARIABLES (TIME SPENT) */
FILE* timing_file;
double time_LCM_init = 0;
double time_initialisation = 0;
double time_threshold_correction = 0;
double time_termination = 0;
double t_init,t_end;
double tic,toc;

// Measure running time
double measureTime(){
  struct rusage t;
  struct timeval tv,ts;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  ts = t.ru_stime;
  return tv.tv_sec + ts.tv_sec + ((double)tv.tv_usec + (double)ts.tv_usec) * 1e-6;
}

// Measure peak memory usage
size_t measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return (size_t)t.ru_maxrss;
}

// Display execution time and memory consumption
void profileCode(){
	size_t peak_memory;

	fprintf(timing_file,"CODE PROFILING\n");
	fprintf(timing_file,"Total execution time: %f (s).\n",t_end-t_init);
	fprintf(timing_file,"\t Time to initialise LCM: %f (s).\n",time_LCM_init);
	fprintf(timing_file,"\t Time to initialise caches: %f (s).\n",time_initialisation);
	fprintf(timing_file,"\t Time to enumerate significant patterns: %f (s).\n",time_threshold_correction);
	fprintf(timing_file,"\t Time to terminate algorithm: %f (s).\n",time_termination);

	peak_memory = measurePeakMemory();
	fprintf(timing_file,"Peak memory consumption: %lld (KB in Linux, B in Mac OS X).\n",peak_memory);

	printf("Total execution time: %f (s).\n",t_end-t_init);
	printf("Peak memory consumption: %lld (KB in Linux, B in Mac OS X).\n",peak_memory);

	// Close timing file
	fclose(timing_file);
}

#endif
