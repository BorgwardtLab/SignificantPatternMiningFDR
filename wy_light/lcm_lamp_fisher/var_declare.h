

/* This file declares variables which will be used across different C files */

// CODE DEPENDENCIES WITH LCM (typedef statements and such)
#include"lib_e.c"
// CODE DEPENDENCIES WITH NEW CODE (typedef statements and such)
#include"transaction_keeping.h"

// VARIABLES DEFINED IN lcm_var.c
extern int **LCM_Ot, **LCM_Os;
extern int LCM_th;
extern ARY LCM_Trsact;
// VARIABLES DEFINED IN transaction_keeping.c
extern BM_TRANS_LIST current_trans;
// VARIABLES DEFINED IN lamp.c
extern int N;
extern int n;
extern long long effective_total_dataset_frq;


