
#include "BFGS_Optimization.h"


int main( int argc, char *argv[] ){

    if(argc != 2){

        std::cout << "./LinearOpticalSimulation [epsilon]" << std::endl << std::endl;

        return 1;

    }

    double epsilon = std::atof( argv[1] );

    BFGS_Optimization optimizer(4e-6,20.0,epsilon);

    for(int i=0;i<100;i++) optimizer.minimize();

    return 0;

}

