#include <vector>
#include "ThresholdOperator.hh"

using namespace std;

TVectorD softthresholdOperator(TectorD *arg, double lambda){
    TVectorD result(arg->GetNrows());
    for (uint i=0;i<result.GetNrows();i++){
        double temp=max(abs((*arg)[i])-lambda,0.0);
        result[i]=sgn((*arg)[i])*temp;
    }
    return result;
}

TVectorD experimentalSoftthresholdOperator(TVectorD *arg, TVectorD *lambdas){
    TVectorD result(arg->GetNrows());
    if (lambdas->GetNrows()==1){
        return softthresholdOperator(args, (*lambdas)[0]);
    }
    for (uint i=0;i<result.GetNrows();i++){
        double temp=max(abs((*arg)[i])-(*lambdas)[i],0.0);
        result[i]=sgn((*arg)[i])*temp;
    }
    return result;
}

class PrimalDualAlgorithm
{
public:
    PrimalDualAlgorithm() = default;
    ~PrimalDualAlgorithm() = default;

    void set_tuning(double, double, double);
    void set_start_values(TVectorD, TVectorD, TVectorD);
    void set_operator(function<TVectorD(TVectorD,double)>);
    TVectorD compute(TVectorD*, TMatrixD*, TMatrixD*, double);

private:
    unsigned int iterations;
    double t1, t2, t3; // tuning parameters
    double gamma; //=1/lambda
    vector<TVectorD> xi, x, x_bar;
    function<TVectorD(TVectorD, double)> threshold_operator;
}