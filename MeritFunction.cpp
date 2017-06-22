#include "MeritFunction.h"

#define PI 3.141592653589793

#define INITIAL_CONDITION_RANDOM_DEGREE 2000




void MeritFunction::setMeritFunction(int intParam){

    diffPhotonNumb = 2;

    std::vector<int> photons,compSubspaceDim;
    std::vector<Eigen::MatrixXi> computationalBasisIn,computationalBasisOut;
    std::vector<Eigen::MatrixXcd> IdealOp;

    photons.resize(diffPhotonNumb);
    compSubspaceDim.resize(diffPhotonNumb);
    IdealOp.resize(diffPhotonNumb);
    computationalBasisIn.resize(diffPhotonNumb);
    computationalBasisOut.resize(diffPhotonNumb);

    photons[0] = 1;
    photons[1] = 2;

    int modes = 6;

    int ancillaPhotons = 1;
    int ancillaModes = 2;

    int measModes = 2;
    int measOutcome = 1;

    compSubspaceDim[0] = 4;
    compSubspaceDim[1] = 8;

    for(int i=0;i<diffPhotonNumb;i++) computationalBasisIn.at(i).resize(compSubspaceDim[i],modes);
    for(int i=0;i<diffPhotonNumb;i++) computationalBasisOut.at(i).resize(compSubspaceDim[i],modes);

    computationalBasisIn[0]   << 0,0,1,0,0,0,
                                 0,0,0,1,0,0,
                                 0,0,0,0,1,0,
                                 0,0,0,0,0,1;

    computationalBasisOut[0]  << 0,0,1,0,0,0,
                                 0,0,0,1,0,0,
                                 0,0,0,0,1,0,
                                 0,0,0,0,0,1;

    computationalBasisIn[1]   << 1,0,1,0,0,0,
                                 1,0,0,1,0,0,
                                 1,0,0,0,1,0,
                                 1,0,0,0,0,1,
                                 0,1,1,0,0,0,
                                 0,1,0,1,0,0,
                                 0,1,0,0,1,0,
                                 0,1,0,0,0,1;

    computationalBasisOut[1]  << 1,0,1,0,0,0,
                                 1,0,0,1,0,0,
                                 1,0,0,0,1,0,
                                 1,0,0,0,0,1,
                                 0,1,1,0,0,0,
                                 0,1,0,1,0,0,
                                 0,1,0,0,1,0,
                                 0,1,0,0,0,1;

    for(int i=0;i<diffPhotonNumb;i++) IdealOp.at(i).resize(compSubspaceDim.at(i),compSubspaceDim.at(i));

    IdealOp[0] << 1.0,0.0,0.0,0.0,
                  0.0,1.0,0.0,0.0,
                  0.0,0.0,1.0,0.0,
                  0.0,0.0,0.0,1.0;

    IdealOp[1] << 0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,
                  1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                  0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,
                  0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,
                  0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,
                  0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,
                  0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,
                  0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0;

    funcDimension = (modes + ancillaModes) * (modes + ancillaModes) + 2*g(ancillaPhotons,ancillaModes);

    La.resize(diffPhotonNumb);

    for(int i=0;i<diffPhotonNumb;i++) La[i].setAncillaAugment(compSubspaceDim[i],ancillaPhotons,ancillaModes);

    Eigen::VectorXd positionTest = Eigen::VectorXd::Random(funcDimension);

    for(int i=0;i<diffPhotonNumb;i++) La[i].setAugmentMatrix(positionTest);

    for(int i=0;i<diffPhotonNumb;i++) La[i].printAugmentMatrix();

    LOCircuit.resize(diffPhotonNumb);

    setLOCircuit(measOutcome,measModes,ancillaPhotons,ancillaModes,computationalBasisIn,computationalBasisOut);

    U.resize(modes + ancillaModes,modes + ancillaModes);

    ASDimension = g(ancillaPhotons,ancillaModes);

    PAULa.resize(diffPhotonNumb);

    return;

}

void MeritFunction::setLOCircuit(int measOutcome,int measModes,int ancillaPhotons,int ancillaModes,std::vector<Eigen::MatrixXi>& compBasisIn,std::vector<Eigen::MatrixXi>& compBasisOut){

    Eigen::MatrixXi ancillaBasis;

    setToFullHilbertSpace(ancillaPhotons,ancillaModes,ancillaBasis);

    Eigen::MatrixXi measBasis;

    setMeasBasis(measOutcome,measModes,measBasis);

    for(int i=0;i<LOCircuit.size();i++){

        Eigen::MatrixXi inBasis,outBasis;

        setInBasis(compBasisIn[i],ancillaBasis,inBasis);

        std::cout << inBasis << std::endl << std::endl;

        setOutBasis(compBasisOut[i],measBasis,outBasis);

        std::cout << outBasis << std::endl << std::endl;

        LOCircuit[i].initializeCircuit(inBasis,outBasis);

    }

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    Eigen::VectorXd a = position.segment(2*ASDimension,U.cols()*U.rows());

    U = genUnitary(a);

    for(int i=0;i<diffPhotonNumb;i++){

        La[i].setAugmentMatrix(position);

        LOCircuit[i].setOmega(U);

        PAULa[i] = LOCircuit[i].omega * La[i].AugmentMatrix;

        std::cout << PAULa[i] << std::endl << std::endl;

        std::cout << (PAULa[i].conjugate().transpose() * PAULa[i]).trace() << std::endl << std::endl;

        std::cout << (PAULa[i] * PAULa[i].conjugate().transpose()).trace() << std::endl << std::endl;

    }

    assert( false );

    return -norm ( (PAULa[0].conjugate().transpose() * PAULa[0]).trace() );

}



void MeritFunction::printReport(Eigen::VectorXd& position){

    return;

}



Eigen::VectorXd MeritFunction::setInitialPosition(){

    U = Eigen::MatrixXcd::Identity(U.rows(),U.cols());

    Eigen::VectorXd output = PI * Eigen::VectorXd::Random(funcDimension);

    std::complex<double> I(0.0,1.0);

    for(int j=0;j<INITIAL_CONDITION_RANDOM_DEGREE;j++){

        Eigen::VectorXd a = Eigen::VectorXd::Random( U.rows() * U.cols() );
        a *= 2000 * PI;
        Eigen::MatrixXcd Utemp( U.rows() , U.cols() );
        Utemp = genUnitary(a);
        U *= Utemp;

    }

    Eigen::MatrixXcd H( U.rows() , U.cols() );

    H = matrixLog(U) / I;

    Eigen::VectorXd a = convertHermittoA(H);

    output.segment( 2*ASDimension , a.size() ) = a;

    return output;

}


void MeritFunction::setOutBasis(Eigen::MatrixXi& compBasisOut,Eigen::MatrixXi& measBasis,Eigen::MatrixXi& outBasis){

    outBasis.resize(compBasisOut.rows(),compBasisOut.cols()+measBasis.cols());

    if(measBasis.cols()>0 ) for(int i=0;i<outBasis.rows();i++) outBasis.block(i,0,1,measBasis.cols()) = measBasis;

    outBasis.block(0,measBasis.cols(),compBasisOut.rows(),compBasisOut.cols()) = compBasisOut;

    return;

}

void MeritFunction::setMeasBasis(int measOutcome,int measModes,Eigen::MatrixXi& measBasis){

    if(measModes == 0){

        measBasis.resize(0,0);
        return;

    }

    int photons = 0;

    int k=0;

    while(true){

        k += g(photons,measModes);

        if(k>measOutcome) break;

        photons++;

    }

    k -= g(photons,measModes);

    int measRow = measOutcome - k;

    setToFullHilbertSpace(photons,measModes,measBasis);

    measBasis = measBasis.block(measRow,0,1,measModes).eval();

    std::ofstream outfile("MeasurementBasis.dat");

    outfile << measBasis << std::endl << std::endl;

    outfile.close();

    return;

}

void MeritFunction::setInBasis(Eigen::MatrixXi& compBasis,Eigen::MatrixXi& ancillaBasis,Eigen::MatrixXi& inBasis){

    if( ancillaBasis.rows() == 0 ){

        inBasis = compBasis;

        return;

    }

    inBasis.resize( compBasis.rows() * ancillaBasis.rows(), compBasis.cols() + ancillaBasis.cols() );

    for(int i=0;i<ancillaBasis.rows();i++){

        inBasis.block(i*compBasis.rows(),ancillaBasis.cols(),compBasis.rows(),compBasis.cols()) = compBasis;

        for(int j=i*compBasis.rows();j<(i+1)*compBasis.rows();j++){

            inBasis.block( j,0,1,ancillaBasis.cols() ) = ancillaBasis.block( i,0,1,ancillaBasis.cols() );

        }

    }

    return;

}


MeritFunction::MeritFunction(){



}


void MeritFunction::setToFullHilbertSpace(const int& subPhotons, const int& subModes,Eigen::MatrixXi& nv){

    if(subPhotons==0 && subModes == 0){

        nv.resize(0,0);

        return;

    }

    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );
    return;;
}

Eigen::MatrixXcd MeritFunction::genUnitary(Eigen::VectorXd a){

    return matrixExp(genHermitian(a));

}

Eigen::MatrixXcd MeritFunction::matrixExp(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();
    std::complex<double> I(0.0,1.0);

                                                                //THIS NEEDS TO BE AS EFFICIENT AS POSSIBLE
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = exp(I*evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result=result+exp(I*evalues(j))*sylvester[j];
    }

    return result;
}


Eigen::VectorXd MeritFunction::convertHermittoA(Eigen::MatrixXcd& H){

    int Hsize = H.rows();

    Eigen::VectorXd output(Hsize*Hsize);

    int rowIndex = 0;

    int outputIndex =0;

    for(int i=0;i<Hsize;i++){

        output(outputIndex) = real(H(rowIndex,rowIndex));
        outputIndex++;

        for(int j=rowIndex+1;j<Hsize;j++){
            output(outputIndex) = sqrt(norm(H(rowIndex,j)));
            outputIndex++;
            output(outputIndex) = arg(H(rowIndex,j));
            outputIndex++;
        }

        rowIndex++;
    }

    return output;

}

Eigen::MatrixXcd MeritFunction::genHermitian(Eigen::VectorXd& a){

    std::complex<double> I(0.0,1.0);
    int Hsize = sqrt(a.size());
    Eigen::MatrixXcd m(Hsize,Hsize);
    int extractIndex=0;                                     //REWRITE THIS FUNCTION IT NEEDS TO BE EFFICIENT- EIGEN SHOULD HAVE A STANDARD ONE

    for(int i=0;i<Hsize;i++){

        m(i,i)=a(extractIndex);
        extractIndex++;

        for(int j=i;j<Hsize;j++){

            if(i!=j){

                m(i,j) = a(extractIndex) * exp(I*a(extractIndex+1));
                m(j,i) = a(extractIndex) * exp(-I*a(extractIndex+1));
                extractIndex++;
                extractIndex++;

            }

        }

    }

    return m;
}


Eigen::MatrixXcd MeritFunction::matrixLog(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXcd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::MatrixXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::MatrixXcd result(matrixSize,matrixSize);
    result = log(evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result = result + log(evalues(j))*sylvester[j];
    }

    return result;
}

