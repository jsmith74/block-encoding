#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED

#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <unsupported/Eigen/KroneckerProduct>
#include "LinearOpticalTransform.h"
#include "AncillaAugment.h"

class MeritFunction{

    public:

        MeritFunction();
        void setMeritFunction(double EPS);
        double f(Eigen::VectorXd& position);
        int funcDimension;
        void printReport(Eigen::VectorXd& position);
        Eigen::VectorXd setInitialPosition();

    private:

        int ASDimension,diffPhotonNumb;

        double globalSuccess,eps;

        std::string filenameMax;

        std::vector<LinearOpticalTransform> LOCircuit;
        std::vector<AncillaAugment> La;
        std::vector<Eigen::MatrixXcd> IdealOp;

        Eigen::MatrixXcd U;

        std::vector<Eigen::MatrixXcd> PAULa;

        std::vector<double> fidelity,successProbability;
        std::vector<int> nonZeroX,nonZeroY;

        void setFidelity(int& i);
        void setSuccessProbability(int& i);

        void setNonZeroXandY();
        void setLOCircuit(int measOutcome,int measModes,int ancillaPhotons,int ancillaModes,std::vector<Eigen::MatrixXi>& compBasisIn,std::vector<Eigen::MatrixXi>& compBasisOut);
        void setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv);
        void setInBasis(Eigen::MatrixXi& compBasis,Eigen::MatrixXi& ancillaBasis,Eigen::MatrixXi& inBasis);
        void setMeasBasis(int measOutcome,int measModes,Eigen::MatrixXi& measBasis);
        void setOutBasis(Eigen::MatrixXi& compBasisOut,Eigen::MatrixXi& measBasis,Eigen::MatrixXi& outBasis);
        void setFullIdealOp(Eigen::MatrixXi& outBasis,Eigen::MatrixXi& compBasisOut,Eigen::MatrixXcd& IdealOp);

        Eigen::MatrixXcd genUnitary(Eigen::VectorXd a);
        Eigen::MatrixXcd matrixLog(Eigen::MatrixXcd X);
        Eigen::MatrixXcd matrixExp(Eigen::MatrixXcd X);
        Eigen::MatrixXcd genHermitian(Eigen::VectorXd& a);
        Eigen::VectorXd convertHermittoA(Eigen::MatrixXcd& H);

        inline int g(const int& n,const int& m);
        inline double doublefactorial(int x);

        void setFilename();
};




inline int MeritFunction::g(const int& n,const int& m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}

inline double MeritFunction::doublefactorial(int x){

    assert(x < 171);

    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}

#endif // MERITFUNCTION_H_INCLUDED
