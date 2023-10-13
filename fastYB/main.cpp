#include <bits/stdc++.h>
#include <getopt.h>

using namespace std;

extern "C"{
void run_fastYB(char* trans_file, char* label_file, double alpha);
void set_seed(int seed);
void set_print_patterns(int flag);
void set_num_perms(int jp);
void set_lambda(double _lambda);
void set_zeta(double zeta);
void no_sstart();
}



int main(int argc, char *argv[] ){
    int seed = -1;
    int k = -1;
    int p = 0;
    char opt;
    while ( ( opt = getopt ( argc, argv, "s:K:l:z:pn" ) ) != -1 ) {
        switch ( opt ) {
            case 's': seed = atoi( optarg ); break;
            case 'K': k = atoi( optarg ); break;
            case 'l': set_lambda(atof( optarg )); break;
            case 'z': set_zeta(atof( optarg )); break;
            case 'p': p = 1; break;
            case 'n': no_sstart(); break;
        }
    }

    if ( argc - optind != 3 ) {
        cout << "Parameters: [-s seed] [-K n_perms] [-l lambda] [-p] transactions labels alpha" << endl;
        return 1;
    }

    double alpha = atof(argv[optind+2]);
    if(p > 0) set_print_patterns(1);
    if(k > 0) set_num_perms(k);
    if(seed > 0) set_seed(seed);
    run_fastYB(argv[optind], argv[optind+1], alpha);
}
