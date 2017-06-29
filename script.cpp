#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h>

int main(){

#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<8;i++){

	usleep(2000000 * omp_get_thread_num());

        std::string commandLine;

        std::stringstream ss;

        double eps = 100 * i * 1e-4;

        ss << eps;

        ss >> commandLine;

        commandLine = "./LinearOpticalSimulation " + commandLine;

        system( commandLine.c_str() );

    }


    return 0;

}
