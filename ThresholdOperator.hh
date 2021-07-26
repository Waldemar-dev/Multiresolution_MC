#ifndef THRESHOLDOPERATOR_HH_
#define THRESHOLDOPERATOR_HH_
#include <vector>
#include "DiscreteWaveletTransformation.hh"

using namespace std;

enum OperatorID
{
    Donoho,
    Munich1,
    Munich2,
    Munich3,
    Positive,
    Test
};

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

class ThresholdOperator
{
public:
    ThresholdOperator(OperatorID, Wavelets);
    ~ThresholdOperator() = default;

    TVectorD apply(TVectorD *x, TVectorD *sigmas);
    TVectorD apply(TVectorD *x, vector<double> *par);
    TVectorD apply(TVectorD*);
    void set_donoho_threshold(double in){lambdaDonoho=in;useGivenLambda=true;}
    double compute_donoho_threshold(TVectorD*);
    double compute_threshold(TVectorD*);
    TVectorD compute_test_threshold(TVectorD*, TVectorD*);

private:
    double median(vector<double>);
    double mad(vector<double>);
    double get_donoho_threshold(vector<double> *wavelet_coefficients, unsigned int n);
    double lambdaDonoho=0;
    bool useGivenLambda=false;
    
    OperatorID id;
    Wavelets waveletID;
};

#endif