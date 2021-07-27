#ifndef CALORIMETERSHOWER_HH_
#define CALORIMETERSHOWER_HH_

#include <sstream>
#include <vector>
#include <map>
#include <random>
#include "ThresholdOperator.hh"
#include "LeastSquaresFit.hh"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/MinimizerOptions.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TVirtualFFT.h"
#include <fstream>
#include <tbb/parallel_for.h>

using namespace std;

class CalorimeterShower
{
public:
    CalorimeterShower(unsigned int,unsigned int,vector<double>*, vector<double>*, Wavelets, OperatorID);
    ~CalorimeterShower() = default;

    void set_x_range(double, double); //sets min and max range
    void set_n_modules(unsigned int);
    void set_n_events(unsigned int in) { n_events = in; }
    void set_par(vector<double> in){par=in;}
    void generate_shower(unsigned int n_events_, unsigned int dim);
    void generate_shower(unsigned int);
    void generate_2D_shower(unsigned int);
    void write_file(){write=true;}
    void histo_to_txt(TH2D*, string);
    void set_fit_attempts(unsigned int in){fit_attempts=in;}
    
    TMatrixD compute_L_dual();
    TMatrixD compute_L();
    TMatrixD compute_L2();
    vector<TMatrixD> get_H_vec(){return H_vec;}
    vector<TMatrixD> get_H_dual_vec(){return H_dual_vec;}
    map<int,double> get_mask(){return mask;}
    map<int,double> get_dual_mask(){return dual_mask;}
    TMatrixD get_high_res_image(){return (*high_res_image);}
    double get_variance(){return variance;}
    double get_tot_error(){return tot_error;}
    double get_tot_error_of_approx(){return tot_error_approx;}

private:
    // general variables
    stringstream file_name;
    bool write;
    unsigned int fit_attempts=100;
    unsigned int dimensions=1;
    //COMPASS and CORAL setup
    unsigned int n_modules;
    unsigned int n_events;
    double module_dimension; //mm
    double x_min;
    double x_max;
    double a1_coral, a2_coral, b1_coral, b2_coral, b3_coral;
    double energy;//GeV
    vector<TH1D> deposition_histograms;
    vector<TH2D> deposition_histograms_2d;
    vector<TMatrixD> low_res_images;
    vector<double> stdev;
    double sigmaE_cell(unsigned int, TVectorD*);
    double sigmaE();

    //Multiresolution variables
    unsigned int enlargement;
    double enlarged_module_dimension;
    double displacement;
    TMatrixD* high_res_image;
    map<unsigned int, vector<double>> targets;
    vector<TMatrixD> H_vec, H_dual_vec;
    map<int, double> mask;
    map<int, double> dual_mask;
    map<int, double> wavelet_mask1, wavelet_mask2, wavelet_mask3;                // = create_wavelet_mask(dual_mask);
    map<int, double> wavelet_dual_mask1, wavelet_dual_mask2, wavelet_dual_mask3; // = create_wavelet_mask(mask);
    vector<map<int, double>> wavelet_masks, wavelet_dual_masks;
    vector<double> par;
    double gamma = 1.0;
    double lambda;

    //Wavelet algorithm
    TVectorD T(TVectorD arg, unsigned int J=2);
    TVectorD algorithm3(TVectorD f, TVectorD g, unsigned int counter, unsigned int counter_max, unsigned int J=2);
    Wavelets wavelet;
    OperatorID thresholdID;
    double tot_error;
    double tot_error_approx;

    //AMP LASSO
    TVectorD amp(TVectorD *g, unsigned int counter, unsigned int counter_max);
    TVectorD amp(TVectorD *g, TVectorD *x, TVectorD *last_x, TVectorD *z, double gamma_threshold, unsigned int counter, unsigned int counter_max);
    double variance;
};
#endif