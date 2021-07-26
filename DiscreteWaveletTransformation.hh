#ifndef DISCRETEWAVELETTRANSFORMATION_HH_
#define DISCRETEWAVELETTRANSFORMATION_HH_

#include <vector>
#include "TH1D.h"

using namespace std;

enum Wavelets
{
    Haar,
    db1,
    db2,
    db3,
    db4,
    db5,
    db6,
    db7,
    db8,
    db9,
    db10,
    db11,
    db12,
    db13,
    db14,
    db15,
    bior1_1,
    bior1_3,
    bior1_5,
    bior2_2,
    bior2_4,
    bior2_6,
    bior2_8,
    bior3_1,
    bior3_3,
    bior3_5,
    bior3_7,
    bior3_9,
    bior4_4,
    bior5_5,
    bior6_8,
    coif1,
    coif2,
    coif3,
    coif4,
    coif5,
    sym2,
    sym3,
    sym4,
    sym5,
    sym6,
    sym7,
    sym8,
    sym9,
    sym10,
    b_spline2,
    haar3,
    wavelet_5_3,
    wavelet_3_5
};

class DiscreteWaveletTransformation
{
public:
    DiscreteWaveletTransformation(Wavelets);
    ~DiscreteWaveletTransformation() = default;

    void analyse(TVectorD *in, vector<double> *out);
    void synthesise(vector<double> *in, TVectorD *out);
    void set_wavelet(Wavelets in) { id = in; }
    void set_J(int in) { J = in; }
    void enable_periodic_extension(bool in){per_extension=in;}
    vector<double> get_flags() { return flags; }
    vector<int> get_length() { return length; }

private:
    void downsamp(vector<double>*, int , vector<double>*);
    void upsamp(vector<double> &sig, int M, vector<double> &sig_u);
    void convfft(TH1D *a, vector<double> &b, vector<double> &c);
    void dwt1(TH1D *in, vector<double> *cA, vector<double> *cD);
    void idwt1(vector<double> &x, vector<double> &cA, vector<double> &cD);
    void filter_coefficients(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
    void vecsum(vector<double> &a, vector<double> &b, vector<double> &c);
    TH1D periodic_extension(TH1D *h_sig, int a, string name);
    vector<double> switch_every_second_sign(vector<double> *in, bool first);
    Wavelets id;
    vector<double> flags;
    vector<int> length;
    int J = 1;
    bool per_extension=true;
};
#endif