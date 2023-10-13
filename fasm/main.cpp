#include <bits/stdc++.h>
#include <getopt.h>

using namespace std;

extern "C"{
void run_by(char* trans_file, char* label_file, double alpha);
void set_seed(int seed);
void set_print_patterns(int flag);
void set_lambda(double _lambda);
void use_stepup();
void no_sstart();
}



int main(int argc, char *argv[] ){
    int p = 0;
    char opt;
    while ( ( opt = getopt ( argc, argv, "punl:s:" ) ) != -1 ) {
        switch ( opt ) {
            case 'p': p = 1; break;
            case 'u': use_stepup(); break;
            case 'n': no_sstart(); break;
            case 'l': set_lambda(atof( optarg )); break;
            case 's': set_seed(atof( optarg )); break;
        }
    }

    if ( argc - optind < 3 || argc - optind > 3 ) {
        cout << "Parameters: [-p] transactions labels alpha" << endl;
        return 1;
    }
  double alpha = atof(argv[optind+2]);
  if(p > 0) set_print_patterns(1);
  run_by(argv[optind], argv[optind+1], alpha);
}
