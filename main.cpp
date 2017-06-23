
#include "BFGS_Optimization.h"


int main( int argc, char *argv[] ){

    BFGS_Optimization optimizer(4e-6,20.0,0);

    for(int i=0;i<1;i++) optimizer.minimize();

    return 0;

}

