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
#include <fstream>
#include <chrono>
#include <math.h>

using namespace std;

enum DetectorType{
    OLGA,
    MAINZ,
    GAMS_ECAL1,
    SHASHLIK,
    GAMSRH,
    GAMS_ECAL2
};

class CalorimeterShower
{
public:
    CalorimeterShower(unsigned int, unsigned int, unsigned int, unsigned int, DetectorType);
    CalorimeterShower(unsigned int,unsigned int,unsigned int, unsigned int, DetectorType, Wavelets, OperatorID);
    ~CalorimeterShower() = default;

    void set_x_range(double, double); //sets min and max range
    void set_y_range(double, double); //sets min and max range
    void set_n_modules_x(unsigned int);
    void set_n_modules_y(unsigned int);
    void set_n_events(unsigned int in) { n_events = in; }
    void set_par(vector<double> in){par=in;}
    void generate_shower(unsigned int n_events_, unsigned int dim);
    void generate_shower(unsigned int);
    void generate_2D_shower(unsigned int);
    void write_file(){write=true;}
    void histo_to_txt(TH2D*, string);
    void set_fit_attempts(unsigned int in){fit_attempts=in;}
    void lednev_fit(string);
    void fix_b1(double);
    void fix_b2(double);
    void read_plot(string, string, string);
    void read_plot(string, string);
    void multiresolution(uint, uint);
    void multiresolution(uint);
    void calibrate_penalty();
    void set_a1_range(double in1, double in2){initial_a1_min=min(in1,in2);initial_a1_max=max(in1,in2);}
    void set_a2_range(double in1, double in2){initial_a2_min=min(in1,in2);initial_a2_max=max(in1,in2);}
    void set_b1_range(double in1, double in2){initial_b1_min=min(in1,in2);initial_b1_max=max(in1,in2);}
    void set_b2_range(double in1, double in2){initial_b2_min=min(in1,in2);initial_b2_max=max(in1,in2);}
    void set_b3_range(double in1, double in2){initial_b3_min=min(in1,in2);initial_b3_max=max(in1,in2);}
    void set_direct_fit(bool in=true){direct_fit=in;}
    void reduce_fit_range_x(unsigned int in){fit_range_x_reduction=in;}
    void reduce_fit_range_y(unsigned int in){fit_range_y_reduction=in;}
    void shift_approximation_in_x(unsigned int in){approximation_x_shift=in;}
    void shift_approximation_in_y(unsigned int in){approximation_y_shift=in;}
    void set_number_of_parameter_pairs(unsigned int in){numbOfArguments=2*in-1;}
    void add_initial_limits_for_a(double in1, double in2);
    void add_initial_limits_for_b(double in1, double in2);
    void compute_wavelet_solution();
    void write_stats_into_file(string);
    void set_real_data_showerparameters_for_shashlik(vector<double>, vector<double>);
    void set_min_lednev_fits(uint in){minLednevFits=in;}
    void limit_variables(){limitedVariables=true;}

    TMatrixD compute_L_dual();
    TMatrixD compute_L();
    // TMatrixD compute_L2(map<uint,vector<uint> > *, map<uint,vector<uint> > *);
    vector<TMatrixD> get_H_vec(){return H_vec;}
    vector<TMatrixD> get_H_dual_vec(){return H_dual_vec;}
    map<int,double> get_mask(){return mask;}
    map<int,double> get_dual_mask(){return dual_mask;}
    TMatrixD get_high_res_image(){return (*high_res_image);}
    double get_variance(){return variance;}
    double get_mse(){return tot_error;}
    double get_mse_of_ncdfs(){return tot_error_ncdf;}
    vector<unsigned int> convert_index_to_bins(unsigned int, unsigned int, unsigned int);
    unsigned int convert_bins_to_index(unsigned int, unsigned int, unsigned int);

private:
    void cut_down_TH2(TH2D*);
    void fill_real_data_parameters();
    TMatrixD Kronecker_product(TMatrixD*, TMatrixD*, map<uint,vector<uint> > *, map<uint,vector<uint> > *);
    void create_x_projection(TH2D*, TH1D*);
    // general variables
    stringstream file_name;
    bool write;
    unsigned int fit_attempts=10;
    unsigned int dimensions=2;
    string toolbox_file_name;
    string toolbox_directory_name="";
    string toolbox_th2d_name;
    bool use_toolbox=false;
    unsigned int fit_range_x_reduction=0;
    unsigned int fit_range_y_reduction=0;
    TVectorD *high_res_approximation=0;
    TH1D *hist_fine_1D=0;
    TH2D *hist_fine_2D=0;
    TH1D *cdf_fine_1D = 0;
    TH2D *cdf_fine_2D = 0;
    bool mc_sim=false;
    double amp_time=0;
    double out_best_chi_squared=DBL_MAX;
    double out_mse_showers=0;
    double out_mse_ncdfs=0;
    DetectorType detector_type;
    map<DetectorType,vector<double> > real_data_a, real_data_b, a_coral, b_coral;
    vector<double> mc_truth_a, mc_truth_b;
    //COMPASS and CORAL setup
    unsigned int n_modules_x, n_modules_y;
    unsigned int n_events;
    double module_dimension; //mm
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double energy = 40.0;//GeV
    vector<TH1D> deposition_histograms;
    vector<TH2D> deposition_histograms_2d;
    vector<TMatrixD> low_res_images;
    vector<double> real_data_a_shashlik={0.885, -0.14};
    vector<double> real_data_b_shashlik={8.104, 55.86, 1.52};
    vector<double> coral_a, coral_b;
    TH1D* real_data_ncdf_x_projection=0;
    double sigmaE_cell(unsigned int, TVectorD*);
    double sigmaE_cell(double);
    double sigmaE();
    void construct_real_data_plots(TFile*);
    string construct_atan_sum(int, vector<string>*);
    string construct_shower_sum(uint, vector<string>*);
    TF2* construct_ncdf_function(vector<double>*, vector<double>*, string);
    TF2* construct_shower_function(vector<double>*, vector<double>*, string);
    void convert_TF2_NCDF_to_TH2D(TF2*, TH2D*);

    //Multiresolution variables
    unsigned int enlargement;
    double enlarged_module_dimension;
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
    double lambda=6.61918e-08;
    unsigned int approximation_x_shift=0;
    unsigned int approximation_y_shift=0;
    map<uint,vector<uint> > out_pairs_1, out_pairs_2, dual_out_pairs_1, dual_out_pairs_2;
    TMatrixD *L=0;
    TMatrixD *L_dual=0;
    TVectorD *epsilon=0;

    //Wavelet algorithm
    TVectorD T(TVectorD arg, unsigned int J=2);
    TVectorD algorithm3(TVectorD f, TVectorD g, unsigned int counter, unsigned int counter_max, unsigned int J=2);
    Wavelets wavelet;
    OperatorID thresholdID;
    double tot_error_ncdf = 0;
    double tot_error = 0;

    //AMP LASSO
    uint fit_attempts_lambda=1000;
    TVectorD amp(TVectorD *g, unsigned int counter, unsigned int counter_max);
    TVectorD amp(TVectorD *g, TVectorD *x, TVectorD *last_x, TVectorD *z, double gamma_threshold, unsigned int counter, unsigned int counter_max);
    double variance;
    unsigned int amp_iterations=3;
    double debias_tolerance_increase=10;
    double debias_tolerance=pow(10,1);

    // Lednev fit
    vector<double> lednev_a, lednev_a_errors, lednev_b, lednev_b_errors;
    double lednev_red_chi_sq, lednev_fit_status, lednev_real_data_mse;
    void compute_2D_NCDF(TH2D*, TH2D*);
    void compute_NCDF(TH1D*, TH1D*);
    void NCDF_difference(vector<double>*, vector<double>*, vector<double>*, vector<double>*, TGraph2DErrors*);
    void NCDF_quotient(vector<double>*, vector<double>*, vector<double>*, vector<double>*, TGraph2DErrors*);
    double mse_ncdf=0;
    double fixed_b1=0;
    double fixed_b2=0;
    unsigned int numbOfArguments = 5;
    unsigned int lednev_fit_attempts=100;
    double lednev_tolerance=pow(10,3);
    bool fix_b1_var=false;
    bool fix_b2_var=false;
    unsigned int correlation_hist_binning=500;
    vector<vector<double> > initial_a, initial_b;
    double initial_a1_min=-4, initial_a1_max=4;
    double initial_b1_min=0.1, initial_b1_max=12;
    double initial_a2_min=-2, initial_a2_max=2; 
    double initial_b2_min=12, initial_b2_max=100;
    double initial_a3_min=-2, initial_a3_max=2;
    double initial_b3_min=100, initial_b3_max=1000;
    double initial_b4_min=1, initial_b4_max=200;
    bool direct_fit=false;
    double tolerance=pow(10,-2);
    uint minLednevFits=100;
    bool limitedVariables=false;
    void sort_lednev_parameters(vector<double>*, vector<double>*,vector<double>*,vector<double>*,TMatrixD *, TH2D*);
    void save_correlations(vector<double>*, vector<double>*,TMatrixD*,TH2D*);
};
#endif