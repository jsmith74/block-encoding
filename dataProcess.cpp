#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

int main(){

    for(int i=0;i<100;i++){

        std::string filename;

        std::stringstream ss;

        double eps = 100 * i * 1e-4;

        ss << eps;

        ss >> filename;

        filename = "Max_Success_" + filename + ".dat";

        std::ifstream infile(filename.c_str());

        if( !infile.is_open() ) continue;

        double epsPrint;

        double fidelity[2];

        double successProbability[2];

        double SPavg;

        infile >> epsPrint;
        infile >> fidelity[0];
        infile >> fidelity[1];
        infile >> successProbability[0];
        infile >> successProbability[1];
        infile >> SPavg;

        infile.close();

        std::ofstream outfile("Max_Success_Total.dat",std::ofstream::app);

        outfile << std::setprecision(16) << epsPrint << "\t" << fidelity[0] << "\t" << fidelity[1] << "\t" << successProbability[0] << "\t" << successProbability[1] << "\t" << SPavg << std::endl;

    }

	return 0;

}
