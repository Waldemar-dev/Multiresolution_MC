#ifndef _AMPLasso_HH_
#define _AMPLasso_HH_

#include <math.h>
#include <vector>
using namespace std;

class AMPLasso {
    public:
    AMPLasso(vector<double>);
    ~AMPLasso()=default;

    private:
    unsigned int n, p;
    double tau;
    double sigma;
    double delta;
    vector<double> beta;
    vector<double> gamma_t;
    vector<double> next_beta(TMatrixD*, TVectorD*);
    vector<double> r,last_r;
    double next_r(TMatrixD*);
    double M(double);
    vector<double> tau;
    void next_tau();
};

#endif