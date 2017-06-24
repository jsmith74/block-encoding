#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <stdlib.h>

int main(){

#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<100;i++){

        std::string commandLine;

        std::stringstream ss;

        double eps = 100 * i * 1e-4;

        ss << eps;

        ss >> commandLine;

        commandLine = "./LinearOpticalSimulation " + commandLine;

        system( commandLine.c_str() );

    }

}
