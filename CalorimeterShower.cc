#include "CalorimeterShower.hh"

bool is_in_interval(double in, double low, double up)
{
    if (in <= up && in >= low)
    {
        return true;
    }
    return false;
}

double threshold_function(double *x, double *par)
{
    double result;
    double gamma = par[0];
    double sigma = par[1];
    double A = 1.0;
    double C = gamma * gamma + sigma;
    double B = -abs(x[0]);
    double D = -abs(x[0]) * gamma * gamma;
    double a = B / A;
    double b = C / A;
    double c = D / A;
    double p = b - a * a / 3.0;
    double q = 2 * pow(a, 3) / 27.0 - a * b / 3.0 + c;
    double delta = pow(q / 2.0, 2) + pow(p / 3.0, 3);
    if (p == 0)
    {
        result = (pow(a * a * a - 27 * c, 1 / 3.0) - a) / 3.0;
    }
    else
    {
        double Gamma = -q / 2.0 * sqrt(27.0 / abs(pow(p, 3)));
        if (delta <= 0)
        {
            double pi = atan(1) * 4;
            double acosTerm = acos(Gamma);
            if (std::isnan(acosTerm) && Gamma >= 1)
            {
                acosTerm = 0;
            }
            if (par[2] == (double)Munich1)
            {
                double eta0 = (acosTerm) / 3.0;
                double solution1 = (2 * sqrt(a * a - 3 * b) * cos(eta0) - a) / 3.0;
                result = solution1;
            }
            else if (par[2] == (double)Munich2)
            {
                double eta1 = (acosTerm + 2 * pi) / 3.0;
                double solution2 = (2 * sqrt(a * a - 3 * b) * cos(eta1) - a) / 3.0;
                result = solution2;
            }
            else if (par[2] == (double)Munich3)
            {
                double eta2 = (acosTerm + 4 * pi) / 3.0;
                double solution3 = (2 * sqrt(a * a - 3 * b) * cos(eta2) - a) / 3.0;
                result = solution3;
            }
        }
        else if (delta > 0 && p < 0)
        {
            double acoshTerm = acosh(abs(Gamma));
            if (std::isnan(acoshTerm) && Gamma <= 1)
            {
                acoshTerm = 0;
            }
            double eta = acoshTerm / 3.0;
            double thresholded_value = (2 * sqrt(a * a - 3 * b) * sgn(q) * cosh(eta) + a) / 3.0;
            result = -thresholded_value;
        }
        else if (delta > 0 && p > 0)
        {
            double eta = (asinh(Gamma)) / 3.0;
            double thresholded_value = (2 * sqrt(3 * b - a * a) * sinh(eta) - a) / 3.0;
            result = thresholded_value;
        }
        else
        {
            abort();
        }
    }
    if (std::isnan(result))
    {
        cout << "result is nan" << endl;
        abort();
    }
    return result;
}

void CalorimeterShower::fill_real_data_parameters()
{
    vector<double> temp_a_olga = {3.668};
    vector<double> temp_b_olga = {2.085, 2.085};
    vector<double> temp_a_mainz = {0.91};
    vector<double> temp_b_mainz = {3.211, 5 * pow(10, -17)};
    vector<double> temp_a_gams_ecal1 = {0.277};
    vector<double> temp_b_gams_ecal1 = {0.701, 5.119};
    vector<double> temp_a_gams = {1.191};
    vector<double> temp_b_gams = {6.664, 43.777};
    vector<double> temp_a_shashlik = {0.885, -0.14};
    vector<double> temp_b_shashlik = {8.104, 55.86, 1.52};
    real_data_a.insert(pair<DetectorType, vector<double>>(OLGA, temp_a_olga));
    real_data_b.insert(pair<DetectorType, vector<double>>(OLGA, temp_b_olga));
    real_data_a.insert(pair<DetectorType, vector<double>>(MAINZ, temp_a_mainz));
    real_data_b.insert(pair<DetectorType, vector<double>>(MAINZ, temp_b_mainz));
    real_data_a.insert(pair<DetectorType, vector<double>>(GAMS_ECAL1, temp_a_gams_ecal1));
    real_data_b.insert(pair<DetectorType, vector<double>>(GAMS_ECAL1, temp_b_gams_ecal1));
    real_data_a.insert(pair<DetectorType, vector<double>>(GAMS_ECAL2, temp_a_gams));
    real_data_b.insert(pair<DetectorType, vector<double>>(GAMS_ECAL2, temp_b_gams));
    real_data_a.insert(pair<DetectorType, vector<double>>(GAMSRH, temp_a_gams));
    real_data_b.insert(pair<DetectorType, vector<double>>(GAMSRH, temp_b_gams));
    real_data_a.insert(pair<DetectorType, vector<double>>(SHASHLIK, temp_a_shashlik));
    real_data_b.insert(pair<DetectorType, vector<double>>(SHASHLIK, temp_b_shashlik));
}

CalorimeterShower::CalorimeterShower(unsigned int in_enlargement, unsigned int in_n_modules_x, unsigned int in_n_modules_y, unsigned int in_dimensions, DetectorType in_det_type)
{
    ROOT::EnableImplicitMT();
    enlargement = in_enlargement;
    dimensions = in_dimensions;
    module_dimension = 38.3; // mm
    enlarged_module_dimension = module_dimension / enlargement;
    set_n_modules_x(in_n_modules_x);
    set_n_modules_y(in_n_modules_y);
    (*L) = compute_L();
    write = false;
    detector_type = in_det_type;
    fill_real_data_parameters();
}

CalorimeterShower::CalorimeterShower(unsigned int in_enlargement, unsigned int in_n_modules_x, unsigned int in_n_modules_y, unsigned int dimensions_, DetectorType in_det_type, Wavelets name, OperatorID id)
{
    ROOT::EnableImplicitMT();
    detector_type = in_det_type;
    fill_real_data_parameters();
    wavelet = name;
    thresholdID = id;
    enlargement = in_enlargement;
    dimensions = dimensions_;
    module_dimension = 38.3; // mm
    enlarged_module_dimension = module_dimension / enlargement;
    set_n_modules_x(in_n_modules_x);
    set_n_modules_y(in_n_modules_y);
    (*L) = compute_L();
    if (enlargement == 2 || enlargement == 4)
    {
        (*L_dual) = compute_L_dual();
    }
    write = false;

    uint n_entries = n_modules_x * enlargement;
    if (in_n_modules_y > 1)
    {
        n_entries *= (n_modules_y * enlargement);
    }
    TMatrixD H(n_entries, n_entries);
    TMatrixD H_dual(H.GetNrows(), H.GetNcols());
    unsigned int index = module_dimension / enlarged_module_dimension;
    for (uint i = 0; i <= index / 2; i++)
    {
        mask[i - 1] = 1.0 / enlargement;
    }
    if (index / 2.0 == (int)(index / 2.0))
    {
        mask[(unsigned int)(index / 2.0) - 1] = 1.0 / 2.0 / enlargement;
    }
    if (enlargement == 2)
    {
        mask[-1] = 1 / 4.0;
        mask[0] = 1 / 2.0;
        mask[1] = 1 / 4.0;
        dual_mask.insert(pair<int, double>(-2, -0.125));
        dual_mask.insert(pair<int, double>(-1, 0.25));
        dual_mask.insert(pair<int, double>(0, 3 / 4.0));
        dual_mask.insert(pair<int, double>(1, 0.25));
        dual_mask.insert(pair<int, double>(2, -0.125));
        wavelet_mask1.insert(pair<int, double>(-1, 0.125));
        wavelet_mask1.insert(pair<int, double>(0, 0.25));
        wavelet_mask1.insert(pair<int, double>(1, -0.75));
        wavelet_mask1.insert(pair<int, double>(2, 0.25));
        wavelet_mask1.insert(pair<int, double>(3, 0.125));
        wavelet_masks.push_back(wavelet_mask1);
        wavelet_dual_mask1.insert(pair<int, double>(0, 0.25));
        wavelet_dual_mask1.insert(pair<int, double>(-1, -0.5)); // not according to formula
        wavelet_dual_mask1.insert(pair<int, double>(-2, 0.25)); // not according to formula
        wavelet_dual_masks.push_back(wavelet_dual_mask1);
    }
    else if (enlargement == 4)
    {
        mask[-2] = 1 / 8.0;
        mask[-1] = 1 / 4.0;
        mask[0] = 1 / 4.0;
        mask[1] = 1 / 4.0;
        mask[2] = 1 / 8.0;
        dual_mask.insert(pair<int, double>(-3, -1 / 16.0));
        dual_mask.insert(pair<int, double>(-2, 1 / 8.0));
        dual_mask.insert(pair<int, double>(-1, 5 / 16.0));
        dual_mask.insert(pair<int, double>(0, 1 / 4.0));
        dual_mask.insert(pair<int, double>(1, 5 / 16.0));
        dual_mask.insert(pair<int, double>(2, 1 / 8.0));
        dual_mask.insert(pair<int, double>(3, -1 / 16.0));
        /*wavelet_mask1.insert(pair<int,double>(-2,1/16.0));
    wavelet_mask1.insert(pair<int,double>(-1,1/8.0));
    wavelet_mask1.insert(pair<int,double>(0,-3/8.0));
    wavelet_mask1.insert(pair<int,double>(1,1/8.0));
    wavelet_mask1.insert(pair<int,double>(2,1/16.0));
    wavelet_mask2.insert(pair<int,double>(-2,5/32.0));
    wavelet_mask2.insert(pair<int,double>(-1,5/16.0));
    wavelet_mask2.insert(pair<int,double>(0,-3/16.0));
    wavelet_mask2.insert(pair<int,double>(1,-3/16.0));
    wavelet_mask2.insert(pair<int,double>(2,-3/32.0));
    wavelet_mask3.insert(pair<int,double>(-6,1/64.0));
    wavelet_mask3.insert(pair<int,double>(-5,1/32.0));
    wavelet_mask3.insert(pair<int,double>(-4,1/32.0));
    wavelet_mask3.insert(pair<int,double>(-3,1/32.0));
    wavelet_mask3.insert(pair<int,double>(-2,-4/64.0));
    wavelet_mask3.insert(pair<int,double>(-1,11/32.0));
    wavelet_mask3.insert(pair<int,double>(0,-5/32.0));
    wavelet_mask3.insert(pair<int,double>(1,-5/32.0));
    wavelet_mask3.insert(pair<int,double>(2,-5/64.0));
    wavelet_dual_mask1.insert(pair<int, double>(0, -2.0));
    wavelet_dual_mask1.insert(pair<int, double>(1, 2.0));
    wavelet_dual_mask2.insert(pair<int, double>(-3, -2));
    wavelet_dual_mask2.insert(pair<int, double>(-2, 4.0));
    wavelet_dual_mask2.insert(pair<int, double>(1, -2.0));
    wavelet_dual_mask3.insert(pair<int, double>(-3, 2.0));
    wavelet_dual_mask3.insert(pair<int, double>(-2, -4.0));
    wavelet_dual_mask3.insert(pair<int, double>(-1, 2.0));*/
        wavelet_mask1.insert(pair<int, double>(-2, -1.0 / 8.0));
        wavelet_mask1.insert(pair<int, double>(-1, -1.0 / 4.0));
        wavelet_mask1.insert(pair<int, double>(0, 0));
        wavelet_mask1.insert(pair<int, double>(1, 1.0 / 4.0));
        wavelet_mask1.insert(pair<int, double>(2, 1.0 / 8.0));
        wavelet_mask2.insert(pair<int, double>(-2, -1.0 / 16.0));
        wavelet_mask2.insert(pair<int, double>(-1, -1.0 / 8.0));
        wavelet_mask2.insert(pair<int, double>(0, 5.0 / 16.0));
        wavelet_mask2.insert(pair<int, double>(1, -1.0 / 4.0));
        wavelet_mask2.insert(pair<int, double>(2, 5.0 / 16.0));
        wavelet_mask2.insert(pair<int, double>(3, -1.0 / 8.0));
        wavelet_mask2.insert(pair<int, double>(4, -1.0 / 16.0));
        wavelet_mask3.insert(pair<int, double>(-2, 1.0 / 16.0));
        wavelet_mask3.insert(pair<int, double>(-1, 1.0 / 8.0));
        wavelet_mask3.insert(pair<int, double>(0, -7.0 / 16.0));
        wavelet_mask3.insert(pair<int, double>(1, 0.0));
        wavelet_mask3.insert(pair<int, double>(2, 7.0 / 16.0));
        wavelet_mask3.insert(pair<int, double>(3, -1.0 / 8.0));
        wavelet_mask3.insert(pair<int, double>(4, -1.0 / 16.0));
        wavelet_masks.push_back(wavelet_mask1);
        wavelet_masks.push_back(wavelet_mask2);
        wavelet_masks.push_back(wavelet_mask3);
        wavelet_dual_mask1.insert(pair<int, double>(-3, -1.0 / 16.0));
        wavelet_dual_mask1.insert(pair<int, double>(-2, 1.0 / 8.0));
        wavelet_dual_mask1.insert(pair<int, double>(-1, 7.0 / 16.0));
        wavelet_dual_mask1.insert(pair<int, double>(0, 0.0));
        wavelet_dual_mask1.insert(pair<int, double>(1, -7.0 / 16.0));
        wavelet_dual_mask1.insert(pair<int, double>(2, -1.0 / 8.0));
        wavelet_dual_mask1.insert(pair<int, double>(3, -1.0 / 16.0));

        wavelet_dual_mask2.insert(pair<int, double>(-1, 1.0 / 8.0));
        wavelet_dual_mask2.insert(pair<int, double>(0, -1.0 / 4.0));
        wavelet_dual_mask2.insert(pair<int, double>(1, 1.0 / 4.0));
        wavelet_dual_mask2.insert(pair<int, double>(2, -1.0 / 4.0));
        wavelet_dual_mask2.insert(pair<int, double>(3, 1.0 / 8.0));

        wavelet_dual_mask3.insert(pair<int, double>(-1, 1.0 / 8.0));
        wavelet_dual_mask3.insert(pair<int, double>(0, -1.0 / 4.0));
        wavelet_dual_mask3.insert(pair<int, double>(1, 0.0));
        wavelet_dual_mask3.insert(pair<int, double>(2, 1.0 / 4.0));
        wavelet_dual_mask3.insert(pair<int, double>(3, -1.0 / 8.0));
        wavelet_dual_masks.push_back(wavelet_dual_mask1);
        wavelet_dual_masks.push_back(wavelet_dual_mask2);
        wavelet_dual_masks.push_back(wavelet_dual_mask3);
    }

    for (uint i = 0; i < wavelet_masks.size(); i++)
    {
        H.Zero();
        H_dual.Zero();
        for (map<int, double>::iterator it = wavelet_masks[i].begin(); it != wavelet_masks[i].end(); it++)
        {
            if (it->first < 0)
            {
                H[0][H.GetNcols() + it->first] = it->second; // needs cross check
            }
            else
            {
                H[0][it->first] = it->second;
            }
        }
        for (map<int, double>::iterator it = wavelet_dual_masks[i].begin(); it != wavelet_dual_masks[i].end(); it++)
        {
            if (it->first < 0)
            {
                H_dual[0][H_dual.GetNcols() + it->first] = it->second;
            }
            else
            {
                H_dual[0][it->first] = it->second;
            }
        }
        uint temp_y_limit = 1;
        if (n_modules_y > 1)
        {
            temp_y_limit = n_modules_y * enlargement;
        }
        for (uint i = 1; i < n_modules_x * enlargement; i++)
        {
            for (uint j = 0; j < temp_y_limit; j++)
            {
                unsigned int temp_index = n_modules_x * enlargement - 1 + j;
                if (temp_index >= n_modules_x * enlargement * temp_y_limit)
                {
                    temp_index = j - 1;
                }
                H[i][j] = H[i - 1][temp_index];
                H_dual[i][j] = H_dual[i - 1][temp_index];
            }
        }
        H_vec.push_back(H);
        H_dual_vec.push_back(H_dual);
    }
}

TMatrixD CalorimeterShower::compute_L_dual()
{
    if (L_dual != 0)
    {
        return (*L_dual);
    }
    TMatrixD result(n_modules_x * enlargement, n_modules_x * enlargement);
    uint temp_y_limit = 1;
    if (n_modules_y > 1)
    {
        temp_y_limit = n_modules_y * enlargement;
    }
    L_dual = new TMatrixD(n_modules_x * enlargement * temp_y_limit, n_modules_x * enlargement * temp_y_limit);
    for (map<int, double>::iterator it = dual_mask.begin(); it != dual_mask.end(); it++)
    {
        if (it->first < 0)
        {
            result[0][n_modules_x * enlargement + it->first] = it->second;
        }
        else
        {
            result[0][it->first] = it->second;
        }
    }
    for (uint i = 1; i < n_modules_x * enlargement; i++)
    {
        for (uint j = 0; j < n_modules_x * enlargement; j++)
        {
            unsigned int temp_index = n_modules_x * enlargement - 1 + j;
            if (temp_index >= n_modules_x * enlargement)
            {
                temp_index = j - 1;
            }
            result[i][j] = result[i - 1][temp_index];
        }
    }
    if (dimensions > 1)
    {
        TMatrixD *new_result = &result;
        for (uint i = 1; i < dimensions; i++)
        {
            (*new_result) = Kronecker_product(new_result, &result, &dual_out_pairs_1, &dual_out_pairs_2);
        }
        (*L_dual) = (*new_result);
        return (*new_result);
    }
    else
    {
        (*L_dual) = result;
        for (uint i = 0; i < L_dual->GetNrows(); i++)
        {
            for (uint j = 0; j < L_dual->GetNcols(); j++)
            {
                if ((*L_dual)[i][j] != 0)
                {
                    if (dual_out_pairs_1.find(i) == dual_out_pairs_1.end())
                    {
                        vector<unsigned int> temp = {j};
                        dual_out_pairs_1.insert(pair<unsigned int, vector<unsigned int>>(i, temp));
                    }
                    else
                    {
                        dual_out_pairs_1[i].push_back(j);
                    }
                    if (dual_out_pairs_2.find(j) == dual_out_pairs_2.end())
                    {
                        vector<unsigned int> temp = {i};
                        dual_out_pairs_2.insert(pair<uint, vector<uint>>(j, temp));
                    }
                    else
                    {
                        dual_out_pairs_2[j].push_back(i);
                    }
                }
            }
        }
    }
    return result;
}

TMatrixD CalorimeterShower::compute_L()
{
    if (L != 0)
    {
        return (*L);
    }
    TMatrixD result(n_modules_x * enlargement, n_modules_x * enlargement);
    uint n_entries = n_modules_x * enlargement;
    if (n_modules_y > 1)
    {
        n_entries *= n_modules_y * enlargement;
    }
    L = new TMatrixD(n_entries, n_entries);
    unsigned int index = module_dimension / enlarged_module_dimension;
    for (uint i = 0; i <= index / 2; i++)
    {
        result[0][i] = 1.0 / enlargement;
        if ((int)(n_modules_x * enlargement - i) < result.GetNrows())
        {
            result[0][n_modules_x * enlargement - i] = 1.0 / enlargement;
        }
    }
    if (index / 2.0 == (int)(index / 2.0))
    {
        result[0][(unsigned int)(index / 2.0)] = 1.0 / 2.0 / enlargement;
        result[0][n_modules_x * enlargement - (unsigned int)(index / 2.0)] = 1.0 / 2.0 / enlargement;
    }
    for (uint i = 1; i < n_modules_x * enlargement; i++)
    {
        for (uint j = 0; j < n_modules_y * enlargement; j++)
        {
            unsigned int temp_index = n_modules_x * enlargement - 1 + j;
            if (temp_index >= n_modules_x * enlargement)
            {
                temp_index = j - 1;
            }
            result[i][j] = result[i - 1][temp_index];
        }
    }
    if (dimensions > 1)
    {
        vector<TMatrixD> temp_results = {result};
        for (uint i = 1; i < dimensions; i++)
        {
            temp_results.push_back(Kronecker_product(&temp_results.back(), &result, &out_pairs_1, &out_pairs_2));
        }
        (*L) = temp_results.back();
        return temp_results.back();
    }
    else
    {
        (*L) = result;
        for (uint i = 0; i < L->GetNrows(); i++)
        {
            for (uint j = 0; j < L->GetNcols(); j++)
            {
                if ((*L)[i][j] != 0)
                {
                    if (out_pairs_1.find(i) == out_pairs_1.end())
                    {
                        vector<unsigned int> temp = {j};
                        out_pairs_1.insert(pair<unsigned int, vector<unsigned int>>(i, temp));
                    }
                    else
                    {
                        out_pairs_1[i].push_back(j);
                    }
                    if (out_pairs_2.find(j) == out_pairs_2.end())
                    {
                        vector<unsigned int> temp = {i};
                        out_pairs_2.insert(pair<uint, vector<uint>>(j, temp));
                    }
                    else
                    {
                        out_pairs_2[j].push_back(i);
                    }
                }
            }
        }
    }
    return result;
}

void CalorimeterShower::set_n_modules_x(unsigned int in)
{
    n_modules_x = in;
    double temp_x_max = n_modules_x * module_dimension / 2.0;
    set_x_range(-temp_x_max, temp_x_max);
}

void CalorimeterShower::set_n_modules_y(unsigned int in)
{
    n_modules_y = in;
    double temp_y_max = n_modules_y * module_dimension / 2.0;
    set_y_range(-temp_y_max, temp_y_max);
}

void CalorimeterShower::set_x_range(double in1, double in2)
{
    x_min = min(in1, in2);
    x_max = max(in1, in2);
}

void CalorimeterShower::set_y_range(double in1, double in2)
{
    y_min = min(in1, in2);
    y_max = max(in1, in2);
}

void CalorimeterShower::generate_shower(unsigned int n_events_)
{
    for (uint i = 0; i < enlargement; i++)
    {
        double x = ((i - (enlargement - 1) / 2.0) * module_dimension / enlargement);
        vector<double> temp2 = {x};
        targets[i] = temp2;
    }
    n_events = n_events_;
    high_res_approximation = new TVectorD(n_modules_x * enlargement * n_modules_y * enlargement);
    TVectorD high_res_control(max(n_modules_x, n_modules_y) * enlargement);
    TFile *my_file = new TFile(file_name.str().c_str(), "recreate");
    TF1 *f1 = new TF1("f1", "([0]*[1]/(pow(x,2)+pow([1],2))+[3]*[2]/(pow(x,2)+pow([2],2))+(1-[0]-[3])*[4]/(pow(x,2)+pow([4],2)))/TMath::Pi()", x_min, x_max);
    f1->SetParameter(0, mc_truth_a[0]); // a1
    f1->SetParameter(1, mc_truth_b[0]); // b1
    f1->SetParameter(2, mc_truth_b[1]); // b2
    f1->SetParameter(3, mc_truth_a[1]); // a2
    f1->SetParameter(4, mc_truth_b[2]); // b3
    f1->SetLineColor(kGreen + 1);
    TF1 *f1Int = new TF1("Int", "([0]*atan((x+[5])/[1]) + [3]*atan((x+[5])/[2]) + (1-[0]-[3])*atan((x+[5])/[4]))/TMath::Pi()+0.5", x_min, x_max);
    f1Int->SetParameter(0, mc_truth_a[0]);          // a1
    f1Int->SetParameter(1, mc_truth_b[0]);          // b1
    f1Int->SetParameter(2, mc_truth_b[1]);          // b2
    f1Int->SetParameter(3, mc_truth_a[1]);          // a2
    f1Int->SetParameter(4, mc_truth_b[2]);          // b3
    f1Int->SetParameter(5, module_dimension / 2.0); // offset
    f1Int->SetLineColor(kGreen + 1);
    if (write)
    {
        f1->Write();
        f1Int->Write();
    }
    deposition_histograms.erase(deposition_histograms.begin(), deposition_histograms.end());
    for (uint i = 0; i < enlargement; i++)
    {
        string name = "Shooting at position " + to_string(i) + " (CORAL);x/mm";
        TH1D hist_temp("Deposition", name.c_str(), max(n_modules_x, n_modules_y), x_min, x_max);
        deposition_histograms.push_back(hist_temp);
    }
    deposition_histograms.shrink_to_fit();
    TH1D *controlHist = new TH1D("control shower", "shower as observed;x/mm", max(n_modules_x, n_modules_y), x_min, x_max);
    controlHist->SetLineColor(kBlue);
    TH1D *controlCDF = new TH1D("control CDF", "NCDF as observed;x/mm", max(n_modules_x, n_modules_y), x_min, x_max);
    controlCDF->SetLineColor(kBlue);
    hist_fine_1D = new TH1D("CORAL Deposition", "High resolution energy deposition (CORAL);x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    hist_fine_1D->SetLineColor(kGreen + 1);
    cdf_fine_1D = new TH1D("NCDF", "Normalized cumulative distribution function (CORAL);x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    cdf_fine_1D->SetLineColor(kGreen + 1);
    TH1D *high_res_hist = new TH1D("g Deposition", "Approximate high resolution energy deposition (g);x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    TH1D *high_res_hist_Lf = new TH1D("L*f Deposition", "Approximate high resolution energy deposition (L*f);x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    for (map<unsigned int, vector<double>>::iterator it = targets.begin(); it != targets.end(); it++)
    {
        for (uint i = 0; i < n_events; i++)
        {
            double r = f1->GetRandom();
            deposition_histograms[it->first].Fill(r + it->second[0]);
            controlHist->Fill(r);
            hist_fine_1D->Fill(r);
        }
        for (int i = 0; i < deposition_histograms[it->first].GetNbinsX(); i++)
        {
            deposition_histograms[it->first].SetBinContent(i + 1, deposition_histograms[it->first].GetBinContent(i + 1) / (double)n_events * energy);
        }
    }
    for (int i = 0; i < controlHist->GetNbinsX(); i++)
    {
        controlHist->SetBinContent(i + 1, controlHist->GetBinContent(i + 1) / ((double)n_events * targets.size()) * energy);
    }
    hist_fine_1D->Scale(energy / (double)hist_fine_1D->Integral());
    TVectorD coral_vec(hist_fine_1D->GetNbinsX());
    for (int i = 0; i < coral_vec.GetNrows(); i++)
    {
        coral_vec[i] = hist_fine_1D->GetBinContent(i + 1);
    }
    TVectorD coral_g((*L) * coral_vec);
    TH1D *coral_g_hist = new TH1D("L*coral", "CORAL cross check;x/mm", coral_g.GetNrows(), x_min, x_max);
    coral_g_hist->SetLineColor(kGreen + 1);
    for (int i = 0; i < coral_g.GetNrows(); i++)
    {
        coral_g_hist->SetBinContent(i + 1, coral_g[i]);
    }
    coral_g_hist->Write();
    if (write)
    {
        controlHist->Write();
    }

    for (int i = 1; i <= controlHist->GetNbinsX(); i++)
    {
        double tempSum = 0;
        for (int j = 1; j <= i; j++)
        {
            tempSum += controlHist->GetBinContent(j);
        }
        controlCDF->SetBinContent(i, tempSum / controlHist->Integral());
    }
    if (write)
    {
        controlCDF->Write();
        for (uint i = 0; i < deposition_histograms.size(); i++)
        {
            deposition_histograms[i].Scale(energy / (double)deposition_histograms[i].Integral());
            for (uint j = 0; j < deposition_histograms[i].GetNbinsX(); j++)
            {
                deposition_histograms[i].SetBinError(j + 1, sigmaE_cell(deposition_histograms[i].GetBinContent(j + 1)));
            }
            deposition_histograms[i].Write();
        }
    }
    cout << endl;

    // rearranging low resolution images into high resolution approximation
    epsilon = new TVectorD(high_res_approximation->GetNrows());
    for (uint i = 0; i < max(n_modules_x, n_modules_y); i++)
    {
        vector<double> bin_entries;
        vector<double> bin_errors;
        for (int j = enlargement - 1; j >= 0; j--)
        {
            bin_entries.push_back(deposition_histograms[j].GetBinContent(i + 1));
            bin_errors.push_back(deposition_histograms[j].GetBinErrorUp(i + 1));
            low_res_images[j][i][0] = deposition_histograms[j].GetBinContent(i + 1);
            low_res_images[j][i][0] = deposition_histograms[j].GetBinErrorUp(i + 1);
        }
        bin_entries.shrink_to_fit();
        bin_errors.shrink_to_fit();
        for (uint j = 0; j < bin_entries.size(); j++)
        {
            high_res_hist->SetBinContent(i * (bin_entries.size()) + j + 1, bin_entries[j]);
            high_res_hist->SetBinError(i * (bin_errors.size()) + j + 1, bin_errors[j]);
            (*high_res_approximation)[i * (bin_entries.size()) + j] = bin_entries[j];
            (*epsilon)[i * (bin_errors.size()) + j] = bin_errors[j];
        }
    }
    TGraphErrors *g_graph = new TGraphErrors(epsilon->GetNrows());
    g_graph->SetTitle("Approximation (g);x/mm;E/GeV");
    double bin_length = (x_max - x_min) / epsilon->GetNrows();
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        double x = (x_max - x_min - bin_length) * i / (double)epsilon->GetNrows() + x_min + bin_length / 2.0;
        g_graph->SetPoint(i, x + bin_length / 2.0, (*high_res_approximation)[i]);
        g_graph->SetPointError(i, bin_length, (*epsilon)[i]);
    }
    g_graph->Write();
    // building NCDF
    for (uint i = 1; i <= max(n_modules_x, n_modules_y) * enlargement; i++)
    {
        double temp = 0;
        for (uint j = 1; j <= i; j++)
        {
            temp += hist_fine_1D->GetBinContent(j);
        }
        cdf_fine_1D->SetBinContent(i, temp / hist_fine_1D->Integral());
    }
    if (write)
    {
        hist_fine_1D->Write();
        high_res_hist->Write();
        cdf_fine_1D->Write();
    }

    TVectorD amp_solution(amp(high_res_approximation, 0, amp_iterations));
    TH1D *amp_histo = new TH1D("amp", "Solution of AMP;x/mm", amp_solution.GetNrows(), x_min, x_max);
    amp_histo->SetLineColor(kRed + 2);
    TH1D *ampIntegral = new TH1D("ampIntegral", "Integral of AMP;x/mm", amp_solution.GetNrows(), x_min, x_max);
    ampIntegral->SetLineColor(kRed + 2);
    double amp_norm = 0;
    for (int i = 0; i < amp_solution.GetNrows(); i++)
    {
        amp_norm += amp_solution[i];
        amp_histo->SetBinContent(i + 1, amp_solution[i]);
    }
    double amp_sum = 0;
    for (int i = 0; i < amp_solution.GetNrows(); i++)
    {
        amp_sum += amp_solution[i];
        ampIntegral->SetBinContent(i + 1, amp_sum / amp_norm);
    }
    amp_histo->Write();
    ampIntegral->Write();
    cout << "variance=" << get_variance() << endl;
    // fit for f
    function<double(const double *)> chisquare_data = chisquare(L, high_res_approximation, epsilon, &out_pairs_1);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(0);
    minimizer->SetTolerance(debias_tolerance);
    // minimizer->SetPrintLevel(1);
    unsigned int numbOfArguments = epsilon->GetNrows();
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments);
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));
    double min_chi_squared = DBL_MAX;
    TH1D *best_f_hist = new TH1D("debiased f", "Debiased f", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    TGraphErrors *best_f_fit = new TGraphErrors(epsilon->GetNrows());
    best_f_fit->SetTitle("debiased f fit;x/mm");
    TH1D *debiased_best_f_integral = new TH1D("debiased_best_f", "debiased f integral;x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    TGraphErrors *min_f_fit = new TGraphErrors(epsilon->GetNrows());
    min_f_fit->SetTitle("debiased f fit;x/mm");
    TH1D *debiased_min_f_integral = new TH1D("debiased_min_f", "debiased f integral;x/mm", max(n_modules_x, n_modules_y) * enlargement, x_min, x_max);
    debiased_min_f_integral->SetLineColor(kBlack);
    TVectorD epsilon_min_chi_sq(max(n_modules_x, n_modules_y) * enlargement);
    TVectorD epsilon_best_chi_sq(max(n_modules_x, n_modules_y) * enlargement);
    double mse_min_f_int = 0;
    out_mse_ncdfs = 0;
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    bool found_fit = false;
    while (!found_fit)
    {
        while (counter < fit_attempts)
        {
            cout << "fit attempt=" << counter + 1 << endl;
            for (int i = 0; i < epsilon->GetNrows(); i++)
            {
                stringstream name;
                name << "f" << to_string(i);
                double start_value = distributions[i](generator);
                if (start_value < 0)
                {
                    start_value = 0;
                }
                minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
                // minimizer->SetVariable(i, name.str().c_str(), start_value, 1e-5);
            }
            minimizer->Minimize();
            counter++;
            if (minimizer->Status() <= 1)
            {
                double bin_length = (x_max - x_min) / epsilon->GetNrows();
                vector<double> args_vec(max(n_modules_x, n_modules_y) * enlargement, 0);
                for (uint i = 0; i < max(n_modules_x, n_modules_y) * enlargement; i++)
                {
                    args_vec[i] = minimizer->X()[i];
                }
                double *args = args_vec.data();
                double temp_chi_sq = chisquare_result(args);
                if (temp_chi_sq < min_chi_squared)
                {
                    double debiased_sum = 0;
                    for (int i = 0; i < epsilon->GetNrows(); i++)
                    {
                        double x = (x_max - x_min - bin_length) * i / (double)epsilon->GetNrows() + x_min + bin_length / 2.0;
                        min_f_fit->SetPoint(i, x + bin_length / 2.0, minimizer->X()[i]);
                        min_f_fit->SetPointError(i, bin_length, minimizer->Errors()[i]);
                        epsilon_min_chi_sq[i] = minimizer->Errors()[i];
                        debiased_sum += minimizer->X()[i];
                        debiased_min_f_integral->SetBinContent(i + 1, debiased_sum);
                    }

                    mse_min_f_int = 0;
                    for (int i = 0; i < debiased_min_f_integral->GetNbinsX(); i++)
                    {
                        debiased_min_f_integral->SetBinContent(i + 1, debiased_min_f_integral->GetBinContent(i + 1) / debiased_sum);
                        mse_min_f_int += pow(cdf_fine_1D->GetBinContent(i + 1) - debiased_min_f_integral->GetBinContent(i + 1), 2);
                    }
                    mse_min_f_int /= epsilon->GetNrows();
                    min_chi_squared = temp_chi_sq;
                }
                if (abs(temp_chi_sq) < abs(out_best_chi_squared))
                {
                    double debiased_sum = 0;
                    for (int i = 0; i < epsilon->GetNrows(); i++)
                    {
                        double x = (x_max - x_min - bin_length) * i / (double)epsilon->GetNrows() + x_min + bin_length / 2.0;
                        best_f_fit->SetPoint(i, x + bin_length / 2.0, minimizer->X()[i]);
                        best_f_fit->SetPointError(i, bin_length, minimizer->Errors()[i]);
                        best_f_hist->SetBinContent(i + 1, minimizer->X()[i]);
                        best_f_hist->SetBinError(i + 1, minimizer->Errors()[i]);
                        epsilon_best_chi_sq[i] = minimizer->Errors()[i];
                        debiased_sum += minimizer->X()[i];
                        debiased_best_f_integral->SetBinContent(i + 1, debiased_sum);
                    }
                    best_f_hist->Scale(energy / (double)best_f_hist->Integral());
                    out_mse_ncdfs = 0;
                    out_mse_showers = 0;
                    for (int i = 0; i < debiased_best_f_integral->GetNbinsX(); i++)
                    {
                        debiased_best_f_integral->SetBinContent(i + 1, debiased_best_f_integral->GetBinContent(i + 1) / debiased_sum);
                        out_mse_ncdfs += pow(cdf_fine_1D->GetBinContent(i + 1) - debiased_best_f_integral->GetBinContent(i + 1), 2);
                        out_mse_showers += pow(best_f_hist->GetBinContent(i + 1) - hist_fine_1D->GetBinContent(i + 1), 2);
                    }
                    out_mse_showers /= debiased_best_f_integral->GetNbinsX();
                    out_mse_ncdfs /= epsilon->GetNrows();
                    out_best_chi_squared = temp_chi_sq;
                    cout << "best_chi_squared=" << out_best_chi_squared << endl;
                }
            }
        }
        if (!found_fit)
        {
            cout << "No fit found. Increasing tolerance by a factor of " << debias_tolerance_increase << "." << endl;
            debias_tolerance *= debias_tolerance_increase;
            minimizer->SetTolerance(debias_tolerance);
            counter = 0;
        }
    }
    best_f_hist->Scale(energy / (double)best_f_hist->Integral());
    best_f_hist->Write("best f hist");
    cout << "lowest red chi^2=" << min_chi_squared << endl;
    min_f_fit->Write("min red chi sq");
    debiased_min_f_integral->Write("min red chi sq int");
    cout << "best red chi^2=" << out_best_chi_squared << endl;
    best_f_fit->Write("best red chi sq");
    debiased_best_f_integral->Write("best red chi sq int");
    cout << "mean squared error of min f int=" << mse_min_f_int << endl;
    cout << "mean squared error of best f int=" << out_mse_ncdfs << endl;
    TGraphErrors *min_red_chi_g = new TGraphErrors();
    TGraphErrors *best_red_chi_g = new TGraphErrors();
    for (int i = 0; i < L->GetNrows(); i++)
    {
        double bin_length = (x_max - x_min) / L->GetNrows();
        double x = (x_max - x_min - bin_length) * i / (double)L->GetNrows() + x_min + bin_length / 2.0;
        double min_chi_sq_value = 0;
        double best_chi_sq_value = 0;
        double min_chi_sq_error = 0;
        double best_chi_sq_error = 0;
        for (int j = 0; j < L->GetNcols(); j++)
        {
            min_chi_sq_value += (*L)[i][j] * min_f_fit->GetY()[j];
            min_chi_sq_error += pow((*L)[i][j] * epsilon_min_chi_sq[j], 2);
            best_chi_sq_value += (*L)[i][j] * best_f_fit->GetY()[j];
            best_chi_sq_error += pow((*L)[i][j] * epsilon_best_chi_sq[j], 2);
        }
        min_chi_sq_error = sqrt(min_chi_sq_error);
        best_chi_sq_error = sqrt(best_chi_sq_error);
        min_red_chi_g->SetPoint(i, x + bin_length / 2.0, min_chi_sq_value);
        min_red_chi_g->SetPointError(i, bin_length / 2.0, min_chi_sq_error);
        best_red_chi_g->SetPoint(i, x + bin_length / 2.0, best_chi_sq_value);
        best_red_chi_g->SetPointError(i, bin_length / 2.0, best_chi_sq_error);
    }
    min_red_chi_g->Write("min_red_chi_sq_g");
    best_red_chi_g->Write("best_red_chi_sq_g");
}

TVectorD CalorimeterShower::T(TVectorD arg, unsigned int J)
{
    TMatrixD L_dual_J((*L_dual));
    TMatrixD L_J((*L));
    if (epsilon->GetNrows() < 1)
    {
        cout << "Errors are not saved." << endl;
        abort();
    }
    for (uint i = 1; i < J; i++)
    {
        L_dual_J *= (*L_dual);
        L_J *= (*L);
    }
    TVectorD result(L_dual_J * L_J * arg);
    ThresholdOperator thresholdOperator(thresholdID, wavelet);
    for (uint j = 0; j < J - 1; j++)
    {
        TMatrixD L_d_j(L_dual->GetNrows(), L_dual->GetNcols());
        L_d_j.UnitMatrix();
        TMatrixD L_j(L->GetNrows(), L->GetNcols());
        L_j.UnitMatrix();
        for (uint i = 0; i < j; i++)
        {
            L_d_j *= (*L_dual);
            L_j *= (*L);
        }
        TVectorD temp(arg.GetNrows());
        temp.Zero();
        if (thresholdID == Donoho)
        {
            lambda = gamma * thresholdOperator.compute_donoho_threshold(&arg);
            thresholdOperator.set_donoho_threshold(lambda);
        }
        for (uint i = 0; i < H_vec.size(); i++)
        {
            TVectorD tempArg(H_vec[i] * L_j * arg);
            if (thresholdID == Donoho || thresholdID == Positive)
            {
                temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg);
            }
            else if (thresholdID == Munich3 || thresholdID == Munich2 || thresholdID == Munich1)
            {
                temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &par);
            }
            else if (thresholdID == Test)
            {
                temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, epsilon);
            }
            else
            {
                cout << "Threshold ID not defined." << endl;
                abort();
            }
        }
        result += L_d_j * temp;
    }
    return result;
}

TVectorD CalorimeterShower::algorithm3(TVectorD f, TVectorD g, unsigned int counter, unsigned int counter_max, unsigned int J)
{
    TVectorD new_f = (*L_dual) * g;
    ThresholdOperator thresholdOperator(thresholdID, wavelet);
    for (uint i = 0; i < H_dual_vec.size(); i++)
    {
        TVectorD temp(H_dual_vec[i] * T(H_vec[i] * f, J));
        new_f += temp;
    }
    if (counter < counter_max)
    {
        new_f = algorithm3(new_f, g, counter + 1, counter_max, J);
    }
    else
    {
        TMatrixD L_dual_J = (*L_dual);
        TMatrixD L_J = (*L);
        for (uint i = 1; i < J; i++)
        {
            L_dual_J *= (*L_dual);
            L_J *= (*L);
        }

        TVectorD result(L_dual_J * L_J * new_f);
        for (uint j = 0; j < J - 1; j++)
        {
            TMatrixD L_d_j(L_dual->GetNrows(), L_dual->GetNcols());
            L_d_j.UnitMatrix();
            TMatrixD L_j(L->GetNrows(), L->GetNcols());
            L_j.UnitMatrix();
            for (uint i = 0; i < j; i++)
            {
                L_d_j *= (*L_dual);
                L_j *= (*L);
            }
            TVectorD temp(new_f.GetNrows());
            temp.Zero();
            if (thresholdID == Donoho)
            {
                lambda = thresholdOperator.compute_donoho_threshold(&new_f);
                thresholdOperator.set_donoho_threshold(lambda);
            }
            for (uint i = 0; i < H_vec.size(); i++)
            {
                TVectorD tempArg(H_vec[i] * L_j * new_f);
                if (thresholdID == Donoho || thresholdID == Positive)
                {
                    temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg);
                }
                else if (thresholdID == Munich1 || thresholdID == Munich2 || thresholdID == Munich3)
                {
                    temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &par);
                }
                else if (thresholdID == Test)
                {
                    temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, epsilon);
                }
                else
                {
                    cout << "Threshold ID not defined." << endl;
                    abort();
                }
            }
            result += L_d_j * temp;
        }
    }
    return new_f;
}

double CalorimeterShower::sigmaE_cell(double in)
{
    double result = sqrt(in / energy) * sigmaE();
    return result;
}

double CalorimeterShower::sigmaE_cell(unsigned int i, TVectorD *in)
{
    double result = sqrt((*in)[i] / energy) * sigmaE();
    return result;
}

double CalorimeterShower::sigmaE()
{
    // double c1 = 0.0114570;
    // double c2 = 0.00145212;
    double c1 = 0.15; // in CORAL: EC02P1__ParamEnergyError
    double c2 = 0.015;
    double c3 = 0.05;

    double result = sqrt(c1 * c1 * energy + c2 * c2 * energy * energy + c3 * c3);
    return result;
}

TVectorD CalorimeterShower::amp(TVectorD *g, TVectorD *x, TVectorD *last_x, TVectorD *z, double gamma_threshold, unsigned int counter, unsigned int counter_max)
{
    if (L->GetNcols() != g->GetNrows())
    {
        abort();
    }
    TMatrixD L_T(L->GetNcols(), L->GetNrows());
    cout << "amp iteration=" << counter << endl;
    L_T.Transpose((*L));
    double next_gamma_threshold = 0;
    double threshold = lambda + gamma_threshold;
    TVectorD next_z((*g) - (*L) * (*x));
    tbb::parallel_for(0, next_z.GetNrows(), [&](int b)
                      {
                          for (int c = 0; c < next_z.GetNrows(); c++)
                          {
                              for (map<uint, vector<uint>>::iterator it = out_pairs_2.begin(); it != out_pairs_2.end(); it++)
                              {
                                  uint j = it->first;
                                  double temp_arg = 0;
                                  for (vector<uint>::iterator d = it->second.begin(); d != it->second.end(); d++)
                                  {
                                      temp_arg += (*L)[*d][j] * (*z)[*d] + (*L)[*d][j] * (*L)[*d][j] * (*last_x)[j];
                                  }
                                  if (temp_arg > gamma_threshold && (*L)[b][j] != 0 && (*L)[c][j] != 0)
                                  {
                                      next_z[b] = next_z[b] + (*L)[b][j] * (*L)[c][j] * (*z)[c];
                                      next_gamma_threshold += (*L)[b][j] * (*L)[c][j];
                                  }
                              }
                          } });
    next_gamma_threshold *= threshold;
    cout << "next_gamma_threshold=" << next_gamma_threshold << endl;
    variance = next_gamma_threshold;
    TVectorD next_x(L_T * next_z);
    for (int i = 0; i < next_x.GetNrows(); i++)
    {
        double temp = 0;
        tbb::parallel_for(0, next_z.GetNrows(), [&](int a)
                          { temp += (*L)[a][i] * next_z[a] + (*L)[a][i] * (*L)[a][i] * (*last_x)[i]; });
        if (temp > variance)
        {
            next_x[i] = temp - variance;
        }
        else
        {
            next_x[i] = 0;
        }
    }
    cout << "next_x=" << endl;
    next_x.Print();

    // logging
    stringstream temp_name;
    temp_name << "AMP_" << to_string(counter + 1);
    TH2D temp_hist(temp_name.str().c_str(), temp_name.str().c_str(), n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    for (uint i = 0; i < next_x.GetNrows(); i++)
    {
        vector<unsigned int> bins = convert_index_to_bins(i, n_modules_x, n_modules_y);
        temp_hist.SetBinContent(bins[0], bins[1], next_x[i]);
        temp_hist.SetBinError(bins[0], bins[1], sqrt(variance));
    }
    amp_results.push_back(temp_hist);
    TVectorD temp_mse_vec((*g) - (*L) * next_x);
    double temp_mse = temp_mse_vec.Norm2Sqr() / (double)temp_mse_vec.GetNrows();
    if (counter == 0)
    {
        best_variance = variance;
        (*best_amp_result) = next_x;
    }
    else if (temp_mse < amp_mse->GetMinimum())
    {
        best_variance = variance;
        (*best_amp_result) = next_x;
    }
    amp_mse->SetPoint(counter, counter + 1, temp_mse);

    if (counter < counter_max)
    {
        TVectorD result(amp(g, &next_x, x, &next_z, next_gamma_threshold, counter + 1, counter_max));
        return result;
    }
    return next_x;
}

TVectorD CalorimeterShower::amp(TVectorD *g, unsigned int counter, unsigned int counter_max) // return amp solution with smalles mean squared error
{
    best_amp_result = new TVectorD(g->GetNrows());
    TVectorD x(g->GetNrows());
    x.Zero();
    TVectorD last_x(x.GetNrows());
    last_x.Zero();
    auto start = std::chrono::high_resolution_clock::now();
    TVectorD result(amp(g, &x, &last_x, g, initial_gamma, counter, counter_max));
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    amp_time = elapsed.count();
    cout << "amp took " << amp_time << " seconds." << endl;
    TVectorD zero_vec(result.GetNrows());
    return (*best_amp_result);
    // if (result == zero_vec)
    // {
    //     cout << "AMP ended in zero vector." << endl;
    //     abort();
    // }
    // return result;
}

void CalorimeterShower::generate_2D_shower(unsigned int n_events_)
{
    n_events = n_events_;
    write = true;
    TFile *my_file = new TFile(file_name.str().c_str(), "recreate");
    TF2 *f2 = construct_shower_function(&mc_truth_a, &mc_truth_b, "MC truth shower");
    TF2 *f2Int = construct_ncdf_function(&mc_truth_a, &mc_truth_b, module_dimension, "MC truth NCDF");

    if (write)
    {
        f2->Write();
        f2Int->Write();
        TH2D *matrix_hist = new TH2D("matrix L", "Matrix L", L->GetNrows(), 0, L->GetNrows(), L->GetNcols(), 0, L->GetNcols());
        for (uint i = 0; i < L->GetNrows(); i++)
        {
            for (uint j = 0; j < L->GetNcols(); j++)
            {
                matrix_hist->SetBinContent(i, j, (*L)[i][j]);
            }
        }
        matrix_hist->Write();
    }

    deposition_histograms_2d.erase(deposition_histograms_2d.begin(), deposition_histograms_2d.end());
    for (uint i = 0; i < enlargement * enlargement; i++)
    {
        string name = "Shooting at position " + to_string(i) + " (CORAL);x/mm";
        stringstream th2d_name;
        th2d_name << "Deposition " << to_string(i);
        TH2D hist_temp(th2d_name.str().c_str(), name.c_str(), n_modules_x, x_min, x_max, n_modules_y, y_min, y_max);
        deposition_histograms_2d.push_back(hist_temp);
    }
    deposition_histograms_2d.shrink_to_fit();
    unsigned int target_counter = 0;
    for (uint i = 0; i < enlargement; i++)
    {
        double x = ((i - (enlargement - 1) / 2.0) * module_dimension / enlargement);
        for (uint j = 0; j < enlargement; j++)
        {
            double y = ((j - (enlargement - 1) / 2.0) * module_dimension / enlargement);
            vector<double> temp2 = {x, y};
            targets[target_counter] = temp2;
            target_counter++;
        }
    }
    TH2D *controlHist = new TH2D("control shower", "shower as observed;x/mm", n_modules_x, x_min, x_max, n_modules_y, y_min, y_max);
    controlHist->SetLineColor(kBlue);
    TH2D *controlCDF = new TH2D("control CDF", "NCDF as observed;x/mm", n_modules_x, x_min, x_max, n_modules_y, y_min, y_max);
    controlCDF->SetLineColor(kBlue);
    hist_fine_2D = new TH2D("CORAL Deposition", "High resolution energy deposition (CORAL);x/mm", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    hist_fine_2D->SetLineColor(kGreen + 1);
    cdf_fine_2D = new TH2D("NCDF", "Normalized cumulative distribution function (CORAL);x/mm", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    cdf_fine_2D->SetLineColor(kGreen + 1);
    TH2D *high_res_hist = new TH2D("g Deposition", "Approximate high resolution energy deposition (g);x/mm", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    TH2D *high_res_hist_Lf = new TH2D("L*f Deposition", "Approximate high resolution energy deposition (L*f);x/mm", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    for (map<unsigned int, vector<double>>::iterator it = targets.begin(); it != targets.end(); it++)
    {
        for (uint i = 0; i < n_events; i++)
        {
            double x, y;
            f2->GetRandom2(x, y);
            deposition_histograms_2d[it->first].Fill(x + it->second[0], y + it->second[1]);
            controlHist->Fill(x, y, 1);
            hist_fine_2D->Fill(x, y, 1);
        }
        for (int i = 0; i < deposition_histograms_2d[it->first].GetNbinsX(); i++)
        {
            for (int j = 0; j < deposition_histograms_2d[it->first].GetNbinsY(); j++)
            {
                deposition_histograms_2d[it->first].SetBinContent(i + 1, j + 1, deposition_histograms_2d[it->first].GetBinContent(i + 1, j + 1) / (double)n_events * energy);
            }
        }
    }
    for (int i = 0; i < controlHist->GetNbinsX(); i++)
    {
        for (int j = 0; j < controlHist->GetNbinsY(); j++)
        {
            controlHist->SetBinContent(i + 1, j + 1, controlHist->GetBinContent(i + 1, j + 1) / ((double)(n_events * targets.size())) * energy);
        }
    }
    hist_fine_2D->Scale(energy / (double)hist_fine_2D->Integral());
    TVectorD coral_vec(hist_fine_2D->GetNbinsX() * hist_fine_2D->GetNbinsY());

    unsigned int coral_vec_counter = 0;
    for (int i = 0; i < hist_fine_2D->GetNbinsX(); i++)
    {
        for (int j = 0; j < hist_fine_2D->GetNbinsY(); j++)
        {
            coral_vec[coral_vec_counter] = hist_fine_2D->GetBinContent(i + 1, j + 1);
            coral_vec_counter++;
        }
    }
    TVectorD coral_g((*L) * coral_vec);
    TH1D *coral_g_hist = new TH1D("L*coral", "CORAL cross check;x/mm", coral_g.GetNrows(), 0, coral_g.GetNrows());
    coral_g_hist->SetLineColor(kGreen + 1);
    for (int i = 0; i < coral_g.GetNrows(); i++)
    {
        coral_g_hist->SetBinContent(i + 1, coral_g[i]);
    }
    coral_g_hist->Write();

    if (write)
    {
        controlHist->Write();
    }
    for (int i = 1; i <= controlHist->GetNbinsX(); i++)
    {
        for (int j = 1; j <= controlHist->GetNbinsY(); j++)
        {
            double tempSum = 0;
            for (int k = 1; k <= i; k++)
            {
                for (int l = 1; l <= j; l++)
                {
                    tempSum += controlHist->GetBinContent(k, l);
                }
            }
            controlCDF->SetBinContent(i, j, tempSum / (double)controlHist->Integral());
        }
    }
    if (write)
    {
        controlCDF->Write();
        for (uint i = 0; i < deposition_histograms_2d.size(); i++)
        {
            deposition_histograms_2d[i].Write();
        }
    }

    // rearranging low resolution images into high resolution approximation (g)
    for (uint k = 0; k < deposition_histograms_2d.size(); k++)
    {
        unsigned int x_offset = floor(k / (double)enlargement);
        unsigned int y_offset = k % enlargement;
        for (uint i = 0; i < n_modules_x; i++)
        {
            double x_bin = enlargement * (1 + i) - x_offset;
            for (uint j = 0; j < n_modules_y; j++)
            {
                double y_bin = enlargement * (j + 1) - y_offset;
                high_res_hist->SetBinContent(x_bin, y_bin, deposition_histograms_2d[k].GetBinContent(i + 1, j + 1) / pow(enlargement, dimensions));
            }
        }
    }
    high_res_hist->Write();
    // rewriting high resolution approximation in vector notation
    high_res_approximation = new TVectorD(L->GetNcols());
    for (int i = 0; i < high_res_hist->GetNbinsX(); i++)
    {
        for (int j = 0; j < high_res_hist->GetNbinsY(); j++)
        {
            (*high_res_approximation)[i * high_res_hist->GetNbinsY() + j] = high_res_hist->GetBinContent(i + 1, j + 1);
        }
    }
    epsilon = new TVectorD(high_res_approximation->GetNrows());
    TGraphErrors *g_graph = new TGraphErrors(epsilon->GetNrows());
    g_graph->SetTitle("Approximation (g);x;E/GeV");
    double bin_length = (x_max - x_min) / (double)(n_modules_x * enlargement);
    unsigned int non_zero_epsilon_entries = 0;
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        (*epsilon)[i] = sigmaE_cell(i, high_res_approximation) / enlargement;
        if ((*epsilon)[i] != 0)
        {
            non_zero_epsilon_entries++;
        }
        g_graph->SetPoint(i, i, (*high_res_approximation)[i]);
        g_graph->SetPointError(i, 0, (*epsilon)[i]);
    }
    g_graph->Write();
    for (uint i = 1; i <= n_modules_x * enlargement; i++)
    {
        for (uint j = 1; j <= n_modules_y * enlargement; j++)
        {
            double temp = 0;
            for (uint k = 1; k <= i; k++)
            {
                for (uint l = 1; l <= j; l++)
                {
                    temp += hist_fine_2D->GetBinContent(k, l);
                }
            }
            cdf_fine_2D->SetBinContent(i, j, temp / hist_fine_2D->Integral());
        }
    }
    if (write)
    {
        hist_fine_2D->Write();
        cdf_fine_2D->Write();
    }
    // Intermediate step: Apply AMP to get LASSO solution and error
    TVectorD amp_solution(amp(high_res_approximation, 0, amp_iterations));
    TH2D *amp_histo = new TH2D("amp", "Solution of AMP", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    amp_histo->SetLineColor(kRed + 2);
    TH2D *ampIntegral = new TH2D("ampIntegral", "Integral of AMP", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    ampIntegral->SetLineColor(kRed + 2);
    for (int i = 0; i < amp_solution.GetNrows(); i++)
    {
        unsigned int x_bin = i / (n_modules_x * enlargement);
        unsigned int y_bin = i % (n_modules_y * enlargement);
        amp_histo->SetBinContent(x_bin + 1, y_bin + 1, amp_solution[i]);
    }
    amp_histo->Write();
    cout << "variance=" << get_variance() << endl;
    // fit for f
    function<double(const double *)> chisquare_data = chisquare(L, high_res_approximation, epsilon, &out_pairs_1);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(0);
    minimizer->SetTolerance(0.1);
    // minimizer->SetPrintLevel(2);
    unsigned int numbOfArguments = epsilon->GetNrows();
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments); // function of type double
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));

    double best_chi_squared = DBL_MAX;
    TGraph2DErrors *best_f_fit = new TGraph2DErrors(epsilon->GetNrows());
    TH2D *temp_hist = new TH2D("best fit histogram", "debiased f fit;x/mm;y/mm;E/GeV", n_modules_x * enlargement, x_min, x_max, n_modules_y * enlargement, y_min, y_max);
    best_f_fit->SetTitle("debiased f fit;x/mm;y/mm;E/GeV");
    TVectorD best_f_vec(n_modules_x * enlargement * n_modules_y * enlargement);
    TVectorD epsilon_best_chi_sq(n_modules_x * enlargement * n_modules_y * enlargement);
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    TH1D *status_hist = new TH1D("status_hist", "Minimizer status", 5, 0, 0);
    while (counter < fit_attempts)
    {
        cout << endl;
        cout << "fit attempt=" << counter << endl;
        for (int i = 0; i < epsilon->GetNrows(); i++)
        {
            stringstream name;
            name << "f" << to_string(i);
            double start_value = distributions[i](generator);
            if (start_value < 0)
            {
                start_value = 0;
            }
            minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
            // minimizer->SetVariable(i, name.str().c_str(), start_value, 1e-5);
        }
        minimizer->Minimize();
        status_hist->Fill(minimizer->Status());
        counter++;
        if (minimizer->Status() <= 1)
        {
            double bin_length = (x_max - x_min) / (n_modules_x * enlargement);
            vector<double> args_vec(epsilon->GetNrows(), 0);
            for (uint i = 0; i < args_vec.size(); i++)
            {
                args_vec[i] = minimizer->X()[i];
            }
            double *args = args_vec.data();
            double temp_chi_sq = chisquare_result(args); // to do
            cout << "chi sq = " << temp_chi_sq << endl;
            if (abs(temp_chi_sq - 1) < abs(best_chi_squared - 1))
            {
                for (int i = 0; i < epsilon->GetNrows(); i++)
                {
                    unsigned int x_bin = i % (n_modules_x * enlargement);
                    unsigned int y_bin = floor(i / (n_modules_y * enlargement));
                    double x = (x_max - x_min - bin_length) * x_bin / (double)(n_modules_x * enlargement) + x_min + bin_length / 2.0;
                    double y = (x_max - x_min - bin_length) * y_bin / (double)(n_modules_y * enlargement) + x_min + bin_length / 2.0;
                    if (isnan(minimizer->X()[i]) || isnan(minimizer->Errors()[i]))
                    {
                        cout << i << "\t" << minimizer->X()[i] << "\t" << minimizer->Errors()[i] << endl;
                    }
                    best_f_fit->SetPoint(i, x + bin_length / 2.0, y + bin_length / 2.0, minimizer->X()[i]);
                    best_f_fit->SetPointError(i, bin_length / 2.0, bin_length / 2.0, minimizer->Errors()[i]);
                    best_f_vec[i] = minimizer->X()[i];
                    epsilon_best_chi_sq[i] = minimizer->Errors()[i];
                    temp_hist->Fill(x, y, minimizer->X()[i]);
                }
                best_chi_squared = temp_chi_sq;
                cout << "best_chi_squared=" << best_chi_squared << endl;
            }
        }
    }

    best_f_fit->GetXaxis()->Set(n_modules_x * enlargement, x_min, x_max);
    best_f_fit->GetYaxis()->Set(n_modules_y * enlargement, x_min, x_max);
    temp_hist->Write("best fit");
    TH2D *debiased_best_f_integral = new TH2D("debiased_best_f", "debiased f NCDF;x/mm;y/mm;", temp_hist->GetNbinsX(), x_min, x_max, temp_hist->GetNbinsY(), x_min, x_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    out_mse_ncdfs = 0;
    for (int i = 1; i <= temp_hist->GetNbinsX(); i++)
    {
        for (int j = 1; j <= temp_hist->GetNbinsY(); j++)
        {
            double temp_sum = 0;
            for (int k = 1; k <= i; k++)
            {
                for (int l = 1; l <= j; l++)
                {
                    temp_sum += temp_hist->GetBinContent(k, l);
                }
            }
            debiased_best_f_integral->SetBinContent(i, j, temp_sum / temp_hist->Integral());
            out_mse_ncdfs += pow(cdf_fine_2D->GetBinContent(i, j) - debiased_best_f_integral->GetBinContent(i, j), 2);
        }
    }
    out_mse_ncdfs /= epsilon->GetNrows();
    compute_2D_NCDF(temp_hist, debiased_best_f_integral);
    best_f_fit->Write("best red chi sq f");
    status_hist->Write();
    debiased_best_f_integral->Write("best red chi sq NCDF");
    cout << "mean squared error of best f int=" << out_mse_ncdfs << endl;
    cout << "red chi squared = " << best_chi_squared << endl;
    TGraphErrors *best_red_chi_g = new TGraphErrors();
    TVectorD best_red_chi_g_vec((*L) * best_f_vec);
    TH1D *deviations_hist = new TH1D("dev", "distribution of nominator with 0 denominator", 100, 0, 0);
    for (int i = 0; i < L->GetNrows(); i++)
    {
        double bin_length = (x_max - x_min) / L->GetNrows();
        double best_chi_sq_value = 0;
        double best_chi_sq_error = 0;
        for (int j = 0; j < L->GetNcols(); j++)
        {
            unsigned int point_number = i * L->GetNcols() + j;
            double x, y, z;
            best_f_fit->GetPoint(point_number, x, y, z);
            best_chi_sq_value += (*L)[i][j] * z;
            best_chi_sq_error += pow((*L)[i][j] * epsilon_best_chi_sq[j], 2);
        }
        best_chi_sq_error = sqrt(best_chi_sq_error);
        best_red_chi_g->SetPoint(i, i, best_red_chi_g_vec[i]);
        best_red_chi_g->SetPointError(i, bin_length / 2.0, best_chi_sq_error);
        if ((*epsilon)[i] == 0)
        {
            deviations_hist->Fill(best_red_chi_g_vec[i]);
        }
    }

    best_red_chi_g->Write("best red chi sq g");
    deviations_hist->Write();
    ofstream stats_file("amp_stats_2d.dat", fstream::app);
    if (stats_file.is_open())
    {
        stats_file << enlargement << "\t" << amp_iterations << "\t" << best_chi_squared << "\t" << out_mse_ncdfs << "\t" << amp_time << endl;
    }
    stats_file.close();
}

void CalorimeterShower::write_stats_into_file(string out_file_name)
{
    ofstream stats_file(out_file_name.c_str(), fstream::app);
    if (stats_file.is_open())
    {
        stats_file << enlargement << "\t" << amp_iterations << "\t" << out_best_chi_squared << "\t" << out_mse_showers << "\t" << out_mse_ncdfs << "\t" << amp_time << endl;
    }
    stats_file.close();
}
// TMatrixD CalorimeterShower::compute_L2(map<uint, vector<uint>> *out_pairs1, map<uint, vector<uint>> *out_pairs2)
// {
//     TMatrixD Lx(compute_L());
//     TMatrixD Ly(Lx);
//     TMatrixD result(Lx.GetNrows() * Ly.GetNrows(), Lx.GetNcols() * Ly.GetNcols());
//     for (int i_x = 0; i_x < Lx.GetNrows(); i_x++)
//     {
//         for (int i_y = 0; i_y < Ly.GetNrows(); i_y++)
//         {
//             for (int j_x = 0; j_x < Lx.GetNcols(); j_x++)
//             {
//                 for (int j_y = 0; j_y < Ly.GetNcols(); j_y++)
//                 {
//                     result[i_x * Ly.GetNrows() + i_y][j_x * Ly.GetNcols() + j_y] = Lx[i_x][j_x] * Ly[i_y][j_y];
//                     if (Lx[i_x][j_x] * Ly[i_y][j_y] != 0)
//                     {
//                         if (out_pairs1->find(i_x * Ly.GetNrows() + i_y) == out_pairs1->end())
//                         {
//                             vector<uint> temp = {(uint)(j_x * Ly.GetNcols() + j_y)};
//                             out_pairs1->insert(pair<uint, vector<uint>>(i_x * Ly.GetNrows() + i_y, temp));
//                         }
//                         else
//                         {
//                             (*out_pairs1)[i_x * Ly.GetNrows() + i_y].push_back(j_x * Ly.GetNcols() + j_y);
//                         }
//                         if (out_pairs2->find(j_x * Ly.GetNcols() + j_y) == out_pairs2->end())
//                         {
//                             vector<uint> temp = {(uint)(i_x * Ly.GetNrows() + i_y)};
//                             out_pairs2->insert(pair<uint, vector<uint>>(j_x * Ly.GetNcols() + j_y, temp));
//                         }
//                         else
//                         {
//                             (*out_pairs2)[j_x * Ly.GetNcols() + j_y].push_back(i_x * Ly.GetNrows() + i_y);
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return result;
// }

TMatrixD CalorimeterShower::Kronecker_product(TMatrixD *in1, TMatrixD *in2, map<uint, vector<uint>> *out_pairs1, map<uint, vector<uint>> *out_pairs2)
{
    TMatrixD result(in1->GetNrows() * in2->GetNrows(), in1->GetNcols() * in2->GetNcols());
    out_pairs1->clear();
    out_pairs2->clear();
    for (int i_x = 0; i_x < in1->GetNrows(); i_x++)
    {
        for (int i_y = 0; i_y < in2->GetNrows(); i_y++)
        {
            for (int j_x = 0; j_x < in1->GetNcols(); j_x++)
            {
                for (int j_y = 0; j_y < in2->GetNcols(); j_y++)
                {
                    result[i_x * in2->GetNrows() + i_y][j_x * in2->GetNcols() + j_y] = (*in1)[i_x][j_x] * (*in2)[i_y][j_y];
                    if ((*in1)[i_x][j_x] * (*in2)[i_y][j_y] != 0)
                    {
                        if (out_pairs1->find(i_x * in2->GetNrows() + i_y) == out_pairs1->end())
                        {
                            vector<uint> temp = {(uint)(j_x * in2->GetNcols() + j_y)};
                            out_pairs1->insert(pair<uint, vector<uint>>(i_x * in2->GetNrows() + i_y, temp));
                        }
                        else
                        {
                            (*out_pairs1)[i_x * in2->GetNrows() + i_y].push_back(j_x * in2->GetNcols() + j_y);
                        }
                        if (out_pairs2->find(j_x * in2->GetNcols() + j_y) == out_pairs2->end())
                        {
                            vector<uint> temp = {(uint)(i_x * in2->GetNrows() + i_y)};
                            out_pairs2->insert(pair<uint, vector<uint>>(j_x * in2->GetNcols() + j_y, temp));
                        }
                        else
                        {
                            (*out_pairs2)[j_x * in2->GetNcols() + j_y].push_back(i_x * in2->GetNrows() + i_y);
                        }
                    }
                }
            }
        }
    }
    return result;
}

void CalorimeterShower::generate_shower(unsigned int n_events_, unsigned int dim)
{
    dimensions = dim;
    mc_sim = true;
    file_name << "mc_tests_x" << to_string(enlargement) << ".root";
    for (uint i = 0; i < enlargement; i++)
    {
        TMatrixD low_res_image(n_modules_x, n_modules_y);
        for (uint j = 0; j < n_modules_x; j++)
        {
            for (uint k = 0; k < n_modules_y; k++)
            {
                low_res_image[j][k] = 0;
            }
        }
        low_res_images.push_back(low_res_image);
    }
    low_res_images.shrink_to_fit();
    if (dim == 1)
    {
        generate_shower(n_events_);
    }
    else if (dim == 2)
    {
        generate_2D_shower(n_events_);
    }
}

void CalorimeterShower::histo_to_txt(TH2D *in, string out_file_name)
{
    ofstream file;
    file.open(out_file_name);
    for (int i = 0; i < in->GetNbinsX(); i++)
    {
        for (int j = 0; j < in->GetNbinsY(); j++)
        {
            file << in->GetBinContent(i + 1, j + 1) << "\t";
        }
        file << endl;
    }
    file.close();
}

void CalorimeterShower::compute_2D_NCDF(TH2D *in, TH2D *out)
{
    if (in->GetNbinsX() != out->GetNbinsX() || in->GetNbinsY() != out->GetNbinsY())
    {
        cout << "Bins do not match" << endl;
        abort();
    }
    for (uint i = 1; i <= in->GetNbinsX(); i++)
    {
        for (uint j = 1; j <= in->GetNbinsY(); j++)
        {
            double temp = 0;
            for (uint k = 1; k <= i; k++)
            {
                for (uint l = 1; l <= j; l++)
                {
                    temp += in->GetBinContent(k, l);
                }
            }
            out->SetBinContent(i, j, temp / in->Integral());
        }
    }
}

void CalorimeterShower::read_plot(string in1, string in2, string in3)
{
    toolbox_file_name = in1;
    toolbox_directory_name = in2;
    toolbox_th2d_name = in3;
    use_toolbox = true;
}

void CalorimeterShower::read_plot(string in1, string in2)
{
    toolbox_file_name = in1;
    toolbox_th2d_name = in2;
    use_toolbox = true;
}

void CalorimeterShower::cut_down_TH2(TH2D *in, TH2D *out)
{
    if (in->GetNbinsX() != out->GetNbinsX() || in->GetNbinsY() != out->GetNbinsY())
    {
        cout << "Bins for cut down do not match." << endl;
        abort();
    }
    for (uint i = 0; i < out->GetNbinsX(); i++)
    {
        for (uint j = 0; j < out->GetNbinsY(); j++)
        {
            if (i < fit_range_x_reduction || i >= in->GetNbinsX() - fit_range_x_reduction)
            {
                out->SetBinContent(i + 1, j + 1, 0);
                out->SetBinError(i + 1, j + 1, 0);
            }
            else if (j < fit_range_y_reduction || j >= in->GetNbinsY() - fit_range_y_reduction)
            {
                out->SetBinContent(i + 1, j + 1, 0);
                out->SetBinError(i + 1, j + 1, 0);
            }
            else
            {
                out->SetBinContent(i + 1, j + 1, in->GetBinContent(i + 1, j + 1));
                out->SetBinError(i + 1, j + 1, in->GetBinErrorUp(i + 1, j + 1));
            }
        }
    }
}

void CalorimeterShower::add_initial_limits_for_a(double in1, double in2)
{
    vector<double> temp = {min(in1, in2), max(in1, in2)};
    initial_a.push_back(temp);
    initial_a.shrink_to_fit();
}

void CalorimeterShower::add_initial_limits_for_b(double in1, double in2)
{
    vector<double> temp = {min(in1, in2), max(in1, in2)};
    initial_b.push_back(temp);
    initial_b.shrink_to_fit();
}

void CalorimeterShower::lednev_fit(string plot_name="best fit;1")
{
    if (numbOfArguments != initial_a.size() + initial_b.size())
    {
        cout << "Not enough initial parameter limits." << endl;
        abort();
    }
    if (use_toolbox)
    {
        file_name.str("");
        if (amp_file != "")
        {
            file_name << get_amp_file();
        }
        else
        {
            file_name << "result_x" << enlargement << ".root";
        }
    }
    TFile *my_file = new TFile(file_name.str().c_str(), "read");
    if (!my_file->IsOpen())
    {
        cout << "File is not open." << endl;
        abort();
    }
    TH2D high_res_hist, high_res_ncdf, cutOut;
    high_res_hist = (*(TH2D *)my_file->Get(plot_name.c_str()));
    high_res_ncdf = (*(TH2D *)my_file->Get("best chi sq NCDF"));
    cutOut = (*(TH2D *)my_file->Get("cut out;1"));
    my_file->Close();
    TH2D *reduced_data = new TH2D("reduced data", "(High resolution) data within defined region;x/mm;y/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    if (direct_fit)
    {
        cut_down_TH2(&cutOut, reduced_data);
    }
    else
    {
        cut_down_TH2(&high_res_hist, reduced_data);
    }
    stringstream new_file_name;
    if (lednev_output_file_name!=""){
        new_file_name<<lednev_output_file_name;
    } else {
        new_file_name << "lednev_fit_x" << enlargement << ".root";
    }
    my_file = new TFile(new_file_name.str().c_str(), "recreate");
    if (direct_fit)
    {
        cutOut.Write("data");
    }
    else
    {
        high_res_hist.Write("data");
    }
    reduced_data->Write("reduced data");
    // build tree
    TTree *lednev_tree = new TTree("Lednev fit results", "Lednev Fit Results");
    TH2D *lednev_shower_hist = new TH2D("Lednev shower", "Lednev shower;x/mm;y/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    TH2D *lednev_ncdf_hist = new TH2D("Lednev NCDF", "Lednev NCDF;x/mm;y/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    TH2D *lednev_correlations = new TH2D("Lednev correlations", "Correlations;a;b", numbOfArguments / 2 + 1, 0.5, numbOfArguments / 2 + 1.5, numbOfArguments / 2 + 1, 0.5, numbOfArguments / 2 + 1.5);
    TH2D *lednev_red_chi_sq_contributions = new TH2D("red chi sq contributions", "Contributions to #chi^{2}_{red.};x/mm;y/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    TH2D *lednev_transform = new TH2D("Transformed Lednev", "Lednev through low-pass filter;x/mm;y/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    TF2 *lednev_ncdf = new TF2();
    TF2 *lednev_shower_function = new TF2();
    TGraph2DErrors *lednev_diff_graph = new TGraph2DErrors();
    TH1D *lednev_ncdf_x_projection = new TH1D("Lednev NCDF x projection", "x projection of the Lednev NCDF;x/mm", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax());
    vector<double> temp_lednev_a(numbOfArguments / 2, 0);
    vector<double> temp_lednev_b(numbOfArguments / 2 + 1, 0);
    vector<double> temp_lednev_a_errors(numbOfArguments / 2, 0);
    vector<double> temp_lednev_b_errors(numbOfArguments / 2 + 1, 0);
    lednev_a = temp_lednev_a;
    lednev_a_errors = temp_lednev_a_errors;
    lednev_b = temp_lednev_b;
    lednev_b_errors = temp_lednev_b_errors;
    degrees_of_freedom = reduced_data->GetNbinsX() * reduced_data->GetNbinsY() - numbOfArguments - 2 * fit_range_x_reduction * reduced_data->GetNbinsY() - 2 * fit_range_y_reduction * reduced_data->GetNbinsX() + 4 * fit_range_x_reduction * fit_range_y_reduction;
    lednev_tree->Branch("reduced_chi_square", &lednev_red_chi_sq);
    lednev_tree->Branch("Fit_Status", &lednev_fit_status);
    lednev_tree->Branch("Shower_histogram", &lednev_shower_hist);
    lednev_tree->Branch("Shower", &lednev_shower_function);
    lednev_tree->Branch("NCDF_histogram", &lednev_ncdf_hist);
    lednev_tree->Branch("NCDF", &lednev_ncdf);
    lednev_tree->Branch("Correlations", lednev_correlations);
    lednev_tree->Branch("Difference_between_Lednev_fit_and_real_data", &lednev_diff_graph);
    lednev_tree->Branch("Mean_Squared_Error_between_Lednev_fit_and_real_data", &lednev_real_data_mse);
    lednev_tree->Branch("Lednev_NCDF_x_projection", &lednev_ncdf_x_projection);
    lednev_tree->Branch("red_chi_sq_contributions", &lednev_red_chi_sq_contributions);
    lednev_tree->Branch("transformed_lednev_solution", &lednev_transform);
    TMatrixD lednev_covariance(numbOfArguments, numbOfArguments);
    for (uint i = 0; i < numbOfArguments - 1; i++)
    {
        stringstream var_name;
        stringstream var_error_name;
        var_error_name << "#Delta ";
        if (i % 2 == 0)
        {
            var_name << "a_{" << to_string(i / 2 + 1) << "}";
            var_error_name << var_name.str();
            lednev_tree->Branch(var_name.str().c_str(), &(lednev_a[i / 2]));
            lednev_tree->Branch(var_error_name.str().c_str(), &lednev_a_errors[i / 2]);
        }
        else
        {
            var_name << "b_{" << to_string(i / 2 + 1) << "}";
            var_error_name << var_name.str();
            lednev_tree->Branch(var_name.str().c_str(), &lednev_b[i / 2]);
            lednev_tree->Branch(var_error_name.str().c_str(), &lednev_b_errors[i / 2]);
        }
    }
    string last_b_name = "b_{" + to_string(numbOfArguments / 2 + 1) + "}";
    lednev_tree->Branch(last_b_name.c_str(), &lednev_b[lednev_b.size() - 1]);
    last_b_name = "#Delta " + last_b_name;
    lednev_tree->Branch(last_b_name.c_str(), &lednev_b_errors[lednev_b_errors.size() - 1]);
    set_x_range(high_res_hist.GetXaxis()->GetXmin() + fit_range_x_reduction * enlarged_module_dimension, high_res_hist.GetXaxis()->GetXmax() - fit_range_x_reduction * enlarged_module_dimension);
    set_y_range(high_res_hist.GetYaxis()->GetXmin() + fit_range_y_reduction * enlarged_module_dimension, high_res_hist.GetYaxis()->GetXmax() - fit_range_y_reduction * enlarged_module_dimension);
    TH2D *data_ncdf = new TH2D("data ncdf", "NCDF of debiased f", high_res_hist.GetNbinsX(), high_res_hist.GetXaxis()->GetXmin(), high_res_hist.GetXaxis()->GetXmax(), high_res_hist.GetNbinsY(), high_res_hist.GetYaxis()->GetXmin(), high_res_hist.GetYaxis()->GetXmax());
    compute_2D_NCDF(&high_res_hist, data_ncdf);
    data_ncdf->Write("data ncdf");
    TH1D *data_ncdf_x_pro = new TH1D("data_ncdf_x_pro", "x projection of the NCDF of data;x/mm", data_ncdf->GetNbinsX(), data_ncdf->GetXaxis()->GetXmin(), data_ncdf->GetXaxis()->GetXmax());
    create_x_projection(data_ncdf, data_ncdf_x_pro);
    data_ncdf_x_pro->Write("data ncdf x projection");
    TF2 *f2 = construct_shower_function(&real_data_a[detector_type], &real_data_b[detector_type], "real data shower");
    f2->Write("real data shower");
    TF2 *f2_int = construct_ncdf_function(&real_data_a[detector_type], &real_data_b[detector_type], module_dimension, "real data NCDF");
    f2_int->Write("real data NCDF");
    function<double(const double *)> chisquare_data = chisquare(reduced_data, energy, numbOfArguments);
    if (direct_fit)
    {
        chisquare_data = chisquare(reduced_data, L, numbOfArguments);
    }
    else if (filterIndices)
    {
        vector<vector<uint>> filtered_indices(filter_indices(radius));
        chisquare_data = chisquare(reduced_data, energy, numbOfArguments, filtered_indices);
    }
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(100000); // 0 stops at ~850!?
    minimizer->SetTolerance(lednev_tolerance);
    // minimizer->SetPrintLevel(2);
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments);
    minimizer->SetFunction(f);
    srand(time(NULL));
    double best_mse_ncdf = DBL_MAX;
    double best_red_chi_sq = DBL_MAX;
    map<double, map<unsigned int, TH1D>> a_dist, b_dist;
    while (lednev_tree->GetEntries() == 0)
    {
        unsigned int counter = 0;
        while (counter < lednev_fit_attempts || (lednev_tree->GetEntries() > 0 && lednev_tree->GetEntries() < minLednevFits))
        {
            cout << endl;
            cout << "fit attempt=" << counter << " (" << min(counter / (double)lednev_fit_attempts, lednev_tree->GetEntries() / (double)minLednevFits) * 100.0 << "\%)" << endl;
            for (uint i = 0; i < numbOfArguments - 1; i++)
            {
                if (i % 2 == 0)
                {
                    double start_value_a = (initial_a[i / 2][1] - initial_a[i / 2][0]) * rand() / (double)RAND_MAX + initial_a[i / 2][0];
                    stringstream var_name;
                    var_name << "a" << to_string(i / 2 + 1);
                    minimizer->SetVariable(i, var_name.str().c_str(), start_value_a, 1e-5);
                }
                else
                {
                    double start_value_b = (initial_b[i / 2][1] - initial_b[i / 2][0]) * rand() / (double)RAND_MAX + initial_b[i / 2][0];
                    stringstream var_name;
                    var_name << "b" << to_string(i / 2 + 1);
                    minimizer->SetLowerLimitedVariable(i, var_name.str().c_str(), start_value_b, 1e-5, 0);
                    // minimizer->SetVariable(i, var_name.str().c_str(), start_value_b, 1e-5);
                }
            }
            double start_value_b = (initial_b[initial_b.size() - 1][1] - initial_b[initial_b.size() - 1][0]) * rand() / (double)RAND_MAX + initial_b[initial_b.size() - 1][0];
            stringstream var_name;
            var_name << "b" << to_string(initial_b.size());
            minimizer->SetLowerLimitedVariable(numbOfArguments - 1, var_name.str().c_str(), start_value_b, 1e-5, 0);
            // minimizer->SetVariable(numbOfArguments - 1, var_name.str().c_str(), start_value_b, 1e-5);
            minimizer->Minimize();
            double temp_edm = minimizer->Edm();
            if (minimizer->Status() == 4)
            {
                while (minimizer->Status() == 4 && minimizer->Edm() >= temp_edm)
                {
                    cout << "Did not reach minimum Edm. Continuing the fit." << endl;
                    minimizer->Minimize();
                }
            }
            counter++;
            if (minimizer->Status() <= 1)
            {
                vector<double> args_vec(numbOfArguments, 0);
                bool skip = false;
                for (uint i = 0; i < numbOfArguments; i++)
                {
                    args_vec[i] = minimizer->X()[i];
                    if (isnan(minimizer->Errors()[i]))
                    {
                        skip = true;
                    }
                    if (i % 2 == 0 && i != numbOfArguments - 1)
                    {
                        temp_lednev_a[i / 2] = minimizer->X()[i];
                        temp_lednev_a_errors[i / 2] = minimizer->Errors()[i];
                    }
                    else
                    {
                        temp_lednev_b[i / 2] = minimizer->X()[i];
                        temp_lednev_b_errors[i / 2] = minimizer->Errors()[i];
                    }
                }
                for (uint i = 0; i < lednev_a.size(); i++)
                {
                    cout << temp_lednev_a[i] << " +/- " << temp_lednev_a_errors[i] << ", " << temp_lednev_b[i] << " +/- " << temp_lednev_b_errors[i] << endl;
                }
                cout << "\t" << temp_lednev_b.back() << " +/- " << temp_lednev_b_errors.back() << endl;
                minimizer->GetCovMatrix(lednev_covariance.GetMatrixArray());
                sort_lednev_parameters(&temp_lednev_a, &temp_lednev_a_errors, &temp_lednev_b, &temp_lednev_b_errors, &lednev_covariance, lednev_correlations);
                lednev_a = temp_lednev_a;
                lednev_a_errors = temp_lednev_a_errors;
                lednev_b = temp_lednev_b;
                lednev_b_errors = temp_lednev_b_errors;
                lednev_fit_status = minimizer->Status();
                // cuts
                if (skip)
                {
                    continue;
                }

                double *args = args_vec.data();
                lednev_red_chi_sq = chisquare_result(args) / (double)(degrees_of_freedom);

                lednev_ncdf = construct_ncdf_function(&lednev_a, &lednev_b, enlarged_module_dimension, "Lednev ncdf");
                lednev_ncdf->SetTitle("Lednev fit NCDF;x/mm;y/mm");
                lednev_ncdf->SetNpx(reduced_data->GetNbinsX());
                lednev_ncdf->SetNpy(reduced_data->GetNbinsY());
                // lednev_ncdf_hist=(TH2D*)lednev_ncdf->CreateHistogram();
                convert_TF2_NCDF_to_TH2D(lednev_ncdf, lednev_ncdf_hist);
                create_x_projection(lednev_ncdf_hist, lednev_ncdf_x_projection);
                lednev_shower_function = construct_shower_function(&lednev_a, &lednev_b, "Lednev shower");
                lednev_shower_function->SetNpx(reduced_data->GetNbinsX());
                lednev_shower_function->SetNpy(reduced_data->GetNbinsY());
                cout << "reduced chi sq=" << lednev_red_chi_sq << endl;
                lednev_shower_function->SetTitle("Lednev fit;x/mm;y/mm");
                lednev_shower_hist = (TH2D *)lednev_shower_function->CreateHistogram();
                lednev_shower_hist->Scale(reduced_data->Integral() / lednev_shower_hist->Integral());
                double temp_max_start[2] = {0, 0}, temp_min_start[2] = {0, 0};
                double temp_max = lednev_ncdf->GetMaximum(temp_max_start);
                double temp_min = lednev_ncdf->GetMinimum(temp_min_start);
                if (abs(temp_max - high_res_ncdf.GetBinContent(high_res_ncdf.GetNbinsX(), high_res_ncdf.GetNbinsY())) > tolerance || abs(temp_min - high_res_ncdf.GetBinContent(1, 1)) > tolerance)
                {
                    cout << "Rejecting solution, because NCDF is between " << temp_min << " and " << temp_max << "." << endl;
                    continue;
                }
                if (a_dist.find(lednev_red_chi_sq) == a_dist.end())
                {
                    map<unsigned int, TH1D> temp_map;
                    stringstream temp_hist_name;
                    stringstream temp_hist_name2;
                    for (uint i = 0; i < initial_a.size() + 1; i++)
                    {
                        temp_hist_name.str("");
                        temp_hist_name2.str("");
                        temp_hist_name << "a" << to_string(i + 1) << "_dist_" << to_string(lednev_red_chi_sq);
                        stringstream temp_hist_name2;
                        temp_hist_name2 << "a_{" << to_string(i + 1) << "} distribution";
                        // TH1D *temp_a_dist = new TH1D(temp_hist_name.str().c_str(), temp_hist_name2.str().c_str(), 100, 0, 0);
                        TH1D temp_a_dist(temp_hist_name.str().c_str(), temp_hist_name2.str().c_str(), 100, 0, 0);
                        temp_map.insert(pair<unsigned int, TH1D>(i, temp_a_dist));
                    }
                    a_dist.insert(pair<double, map<unsigned int, TH1D>>(lednev_red_chi_sq, temp_map));
                    map<unsigned int, TH1D> temp_map2;
                    for (uint i = 0; i < initial_b.size(); i++)
                    {
                        temp_hist_name.str("");
                        temp_hist_name2.str("");
                        temp_hist_name << "b" << to_string(i + 1) << "_dist_" << to_string(lednev_red_chi_sq);
                        temp_hist_name2 << "b_{" << to_string(i + 1) << "} distribution";
                        // TH1D *temp_b_dist = new TH1D(temp_hist_name.str().c_str(), temp_hist_name2.str().c_str(), 100, 0, 0);
                        TH1D temp_b_dist(temp_hist_name.str().c_str(), temp_hist_name2.str().c_str(), 100, 0, 0);
                        temp_map2.insert(pair<unsigned int, TH1D>(i, temp_b_dist));
                    }
                    b_dist.insert(pair<double, map<unsigned int, TH1D>>(lednev_red_chi_sq, temp_map2));
                }

                for (uint i = 0; i < lednev_a.size(); i++)
                {
                    a_dist[lednev_red_chi_sq][i].Fill(lednev_a[i]);
                    b_dist[lednev_red_chi_sq][i].Fill(lednev_b[i]);
                }
                b_dist[lednev_red_chi_sq][lednev_a.size()].Fill(lednev_b.back());
                lednev_real_data_mse = 0;
                NCDF_difference(&lednev_a, &lednev_b, &lednev_a_errors, &lednev_b_errors, lednev_diff_graph);
                lednev_diff_graph->SetTitle("Difference between fit and real data NCDFs;x/mm;y/mm");
                plot_red_chi_sq_contributions_per_bin(reduced_data, lednev_shower_hist, lednev_red_chi_sq_contributions);
                transform_lednev_solution(lednev_shower_hist, lednev_transform);
                lednev_tree->Fill();
            }
            else
            {
                cout << "Failed with status " << minimizer->Status() << "." << endl;
            }
        }
        if (lednev_tree->GetEntries() == 0)
        {
            lednev_tolerance *= 10;
            cout << "Trying Lednev fit with tolerance " << lednev_tolerance << endl;
            minimizer->SetTolerance(lednev_tolerance);
            counter = 0;
        }
    }
    for (map<double, map<unsigned int, TH1D>>::iterator it_a = a_dist.begin(); it_a != a_dist.end(); it_a++)
    {
        for (map<unsigned int, TH1D>::iterator it2 = it_a->second.begin(); it2 != it_a->second.end(); it2++)
        {
            it2->second.Write();
        }
        for (map<unsigned int, TH1D>::iterator it2 = b_dist[it_a->first].begin(); it2 != b_dist[it_a->first].end(); it2++)
        {
            it2->second.Write();
        }
    }
    TTree *sorted_lednev_tree = (TTree *)lednev_tree->CloneTree(0);
    Int_t nb_idx = lednev_tree->BuildIndex("reduced_chi_square");
    TTreeIndex *att_index = (TTreeIndex *)lednev_tree->GetTreeIndex();
    if (att_index == 0)
    {
        lednev_tree->Print();
        abort();
    }
    for (Long64_t i = 0; i < att_index->GetN(); i++)
    {
        lednev_tree->GetEntry(att_index->GetIndex()[i]);
        sorted_lednev_tree->Fill();
    }
    sorted_lednev_tree->Write();
}

void CalorimeterShower::NCDF_difference(vector<double> *a, vector<double> *b, vector<double> *a_error, vector<double> *b_error, TGraph2DErrors *out)
{
    TGraph2DErrors *diff_graph = new TGraph2DErrors();
    unsigned int counter = 0;
    for (uint i = 0; i < n_modules_x * enlargement; i++)
    {
        double x = (x_max - enlarged_module_dimension - x_min) * i / (double)(n_modules_x * enlargement) + x_min + enlarged_module_dimension / 2.0;
        for (uint j = 0; j < n_modules_y * enlargement; j++)
        {
            double y = (y_max - enlarged_module_dimension - y_min) * j / (double)(n_modules_y * enlargement) + y_min + enlarged_module_dimension / 2.0;
            double function_value = 0;
            double last_a = 1;
            double last_a_error = 0;
            double error = 0;
            for (uint k = 0; k < a->size(); k++)
            {
                last_a -= a->at(k);
                function_value += a->at(k) * (atan(x / b->at(k)) + atan(y / b->at(k)) + atan(x * y / b->at(k) / sqrt(b->at(k) * b->at(k) + x * x + y * y)));
                error += pow(a_error->at(k) * (atan(x / b->at(k)) + atan(y / b->at(k)) + atan(x * y / b->at(k) / sqrt(b->at(k) * b->at(k) + x * x + y * y))), 2);
                double term1 = a->at(k) * (x / (b->at(k) * b->at(k) + x * x) * b_error->at(k));
                double term2 = a->at(k) * (y / (b->at(k) * b->at(k) + y * y) * b_error->at(k));
                double term3 = -a->at(k) * (x * y * (2 * b->at(k) * b->at(k) + x * x + y * y) / b->at(k) / b->at(k) / pow(b->at(k) * b->at(k) + x * x + y * y, 1.5) / (1 + pow(x * y / b->at(k) / sqrt(x * x + y * y + b->at(k) * b->at(k)), 2))) * b_error->at(k);
                error += pow(term1 + term2 + term3, 2);
                last_a_error += pow(a_error->at(k), 2);
            }
            last_a_error = sqrt(last_a_error);
            function_value += last_a * (atan(x / b->back()) + atan(y / b->back()) + atan(x * y / b->back() / sqrt(b->back() * b->back() + x * x + y * y)));
            error += pow(last_a_error * (atan(x / b->back()) + atan(y / b->back()) + atan(x * y / b->back() / sqrt(b->back() * b->back() + x * x + y * y))), 2);
            double term1 = last_a * (x / (b->back() * b->back() + x * x) * b_error->back());
            double term2 = last_a * (y / (b->back() * b->back() + y * y) * b_error->back());
            double term3 = -last_a * (x * y * (2 * b->back() * b->back() + x * x + y * y) / b->back() / b->back() / pow(b->back() * b->back() + x * x + y * y, 1.5) / (1 + pow(x * y / b->back() / sqrt(x * x + y * y + b->back() * b->back()), 2))) * b_error->back();
            error += pow(term1 + term2 + term3, 2);
            if (error < 0)
            {
                cout << "error is negative" << endl;
                abort();
            }
            error = sqrt(error);
            double last_a_RD = 1;
            double term = 0;
            for (uint k = 0; k < real_data_a[detector_type].size(); k++)
            {
                last_a_RD -= real_data_a[detector_type][k];
                term += real_data_a[detector_type][k] * (atan(x / real_data_b[detector_type][k]) + atan(y / real_data_b[detector_type][k]) + atan(x * y / real_data_b[detector_type][k] / sqrt(real_data_b[detector_type][k] * real_data_b[detector_type][k] + x * x + y * y)));
            }
            double last_RD_term = last_a_RD * (atan(x / real_data_b[detector_type].back()) + atan(y / real_data_b[detector_type].back()) + atan(x * y / real_data_b[detector_type].back() / sqrt(real_data_b[detector_type].back() * real_data_b[detector_type].back() + x * x + y * y)));
            function_value -= (term + last_RD_term);
            function_value /= (2.0 * TMath::Pi());
            lednev_real_data_mse += pow(function_value, 2);
            diff_graph->SetPoint(counter, x, y, function_value);
            diff_graph->SetPointError(counter, enlarged_module_dimension / 2.0, enlarged_module_dimension / 2.0, error);
            counter++;
        }
    }
    lednev_real_data_mse /= counter;
    (*out) = (*diff_graph);
}

void CalorimeterShower::NCDF_quotient(vector<double> *a, vector<double> *b, vector<double> *a_error, vector<double> *b_error, TGraph2DErrors *out)
{
    TGraph2DErrors *out_graph = new TGraph2DErrors();
    unsigned int counter = 0;
    for (uint i = 0; i < n_modules_x * enlargement; i++)
    {
        double x = (x_max - enlarged_module_dimension - x_min) * i / (double)(n_modules_x * enlargement) + x_min + enlarged_module_dimension / 2.0;
        for (uint j = 0; j < n_modules_y * enlargement; j++)
        {
            double y = (y_max - enlarged_module_dimension - y_min) * j / (double)(n_modules_y * enlargement) + y_min + enlarged_module_dimension / 2.0;
            double function_value = 0;
            double last_a = 1;
            double last_a_error = 0;
            double error = 0;
            for (uint k = 0; k < a->size(); k++)
            {
                last_a -= a->at(k);
                function_value += a->at(k) * (atan(x / b->at(k)) + atan(y / b->at(k)) + atan(x * y / b->at(k) / sqrt(b->at(k) * b->at(k) + x * x + y * y)));
                error += pow(a_error->at(k) * (atan(x / b->at(k)) + atan(y / b->at(k)) + atan(x * y / b->at(k) / sqrt(b->at(k) * b->at(k) + x * x + y * y))), 2);
                double term1 = a->at(k) * (x / (b->at(k) * b->at(k) + x * x) * b_error->at(k));
                double term2 = a->at(k) * (y / (b->at(k) * b->at(k) + y * y) * b_error->at(k));
                double term3 = -a->at(k) * (x * y * (2 * b->at(k) * b->at(k) + x * x + y * y) / b->at(k) / b->at(k) / pow(b->at(k) * b->at(k) + x * x + y * y, 1.5) / (1 + pow(x * y / b->at(k) / sqrt(x * x + y * y + b->at(k) * b->at(k)), 2))) * b_error->at(k);
                error += pow(term1 + term2 + term3, 2);
                last_a_error += pow(a_error->at(k), 2);
                if (isnan(a_error->at(k)) || isnan(b_error->at(k)))
                {
                    return;
                }
            }
            last_a_error = sqrt(last_a_error);
            function_value += last_a * (atan(x / b->back()) + atan(y / b->back()) + atan(x * y / b->back() / sqrt(b->back() * b->back() + x * x + y * y)));
            function_value /= (2.0 * TMath::Pi());
            function_value += 0.25;
            error += pow(last_a_error * (atan(x / b->back()) + atan(y / b->back()) + atan(x * y / b->back() / sqrt(b->back() * b->back() + x * x + y * y))), 2);
            double term1 = last_a * (x / (b->back() * b->back() + x * x) * b_error->back());
            double term2 = last_a * (y / (b->back() * b->back() + y * y) * b_error->back());
            double term3 = -last_a * (x * y * (2 * b->back() * b->back() + x * x + y * y) / b->back() / b->back() / pow(b->back() * b->back() + x * x + y * y, 1.5) / (1 + pow(x * y / b->back() / sqrt(x * x + y * y + b->back() * b->back()), 2))) * b_error->back();
            error += pow(term1 + term2 + term3, 2);
            double last_a_coral = 1;
            double coral_term = 0;
            for (uint k = 0; k < a_coral[detector_type].size(); k++)
            {
                last_a_coral -= a_coral[detector_type][k];
                coral_term += a_coral[detector_type][k] * (atan(x / b_coral[detector_type][k]) + atan(y / b_coral[detector_type][k]) + atan(x * y / b_coral[detector_type][k] / sqrt(b_coral[detector_type][k] * b_coral[detector_type][k] + x * x + y * y)));
            }
            double last_coral_term = last_a_coral * (atan(x / b_coral[detector_type].back()) + atan(y / b_coral[detector_type].back()) + atan(x * y / b_coral[detector_type].back() / sqrt(b_coral[detector_type].back() * b_coral[detector_type].back() + x * x + y * y)));
            function_value /= ((coral_term + last_coral_term) / (2.0 * TMath::Pi()) + 0.25);
            error = error / ((coral_term + last_coral_term) / (2.0 * TMath::Pi()) + 0.25);
            if (isnan(error))
            {
                cout << "error is nan" << endl;
                continue;
            }
            else if (error < 0)
            {
                cout << "error is negative" << endl;
                if (0 > ((coral_term + last_coral_term) / (2.0 * TMath::Pi()) + 0.25))
                {
                    cout << "denominator is negative (" << (coral_term + last_coral_term) / (2.0 * TMath::Pi()) + 0.25 << ")" << endl;
                }
                else
                {
                    cout << "nominator is negative" << endl;
                }
                continue;
            }
            error = sqrt(error);
            out_graph->SetPoint(counter, x, y, function_value);
            out_graph->SetPointError(counter, enlarged_module_dimension / 2.0, enlarged_module_dimension / 2.0, error);
            counter++;
        }
    }
    (*out) = (*out_graph);
}

void CalorimeterShower::fix_b1(double in)
{
    fix_b1_var = true;
    fixed_b1 = in;
}

void CalorimeterShower::fix_b2(double in)
{
    fix_b2_var = true;
    fixed_b2 = in;
}

void CalorimeterShower::multiresolution(uint n_cells_x)
{
    amp_mse = new TGraph(amp_iterations);
    amp_mse->SetTitle(";Iteration;MSE");
    if (use_toolbox)
    {
        cout << "Cannot use a toolbox file for a one dimensional shower profile." << endl;
        abort();
    }
    if (!mc_sim)
    {
        cout << "MC simulation disabled." << endl;
        abort();
    }
    TH1D *chi_sq_hist = new TH1D("chi sq", "chi squared distribution", 500, 0, 0);
    // TH2D *approximation = 0;
    stringstream output_file_name;
    output_file_name << "result_x" << enlargement << ".root";
    TFile *my_file = new TFile(output_file_name.str().c_str(), "recreate");

    TVectorD amp_solution(amp(high_res_approximation, 0, amp_iterations));
    TH1D *amp_histo = new TH1D("amp", "Solution of AMP", (n_cells_x * 2 + 1) * enlargement, x_min, x_max);
    amp_histo->SetLineColor(kRed + 2);
    for (int i = 0; i < amp_solution.GetNrows(); i++)
    {
        amp_histo->SetBinContent(i + 1, amp_solution[i]);
    }
    amp_histo->Write("amp");
    for (uint i = 0; i < amp_results.size(); i++)
    {
        amp_results[i].Write();
    }
    amp_results.clear();
    amp_mse->Write("mse");
    // debiasing
    function<double(const double *)> chisquare_data = chisquare(L, high_res_approximation, epsilon, &out_pairs_1);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(0);
    minimizer->SetTolerance(0.1);
    // minimizer->SetPrintLevel(2);
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, epsilon->GetNrows()); // function of type double
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));
    double best_chi_squared = DBL_MAX;
    TH1D *temp_hist = new TH1D("best fit histogram", "debiased f fit;x/mm;y/mm;E/GeV", (n_cells_x * 2 + 1) * enlargement, x_min, x_max);
    TVectorD best_f_vec((n_cells_x * 2 + 1) * enlargement);
    TVectorD epsilon_best_chi_sq((n_cells_x * 2 + 1) * enlargement);
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    if (get_variance() == 0)
    {
        cout << "0 variance. Abort." << endl;
        abort();
    }
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    TH1D *status_hist = new TH1D("status_hist", "Minimizer status", 5, 0, 0);
    double bin_length = module_dimension / ((double)enlargement);
    while (counter < fit_attempts || best_chi_squared > max_debias_chi_sq)
    {
        cout << endl;
        cout << "fit attempt=" << counter + 1 << endl;
        for (int i = 0; i < epsilon->GetNrows(); i++)
        {
            stringstream name;
            name << "f" << to_string(i);
            double start_value = distributions[i](generator);
            if (start_value < 0)
            {
                start_value = 0;
            }
            minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
            // minimizer->SetVariable(i, name.str().c_str(), start_value, 1e-5);
        }
        minimizer->Minimize();
        double temp_edm = minimizer->Edm();
        if (minimizer->Status() == 4)
        {
            while (temp_edm >= minimizer->Edm() && minimizer->Status() == 4)
            {
                temp_edm = minimizer->Edm();
                minimizer->Minimize();
            }
        }
        status_hist->Fill(minimizer->Status());
        counter++;
        if (minimizer->Status() <= 1)
        {
            vector<double> args_vec(epsilon->GetNrows(), 0);
            for (uint i = 0; i < args_vec.size(); i++)
            {
                args_vec[i] = minimizer->X()[i];
            }
            double *args = args_vec.data();
            double temp_chi_sq = chisquare_result(args);
            cout << "chi sq = " << temp_chi_sq << endl;
            chi_sq_hist->Fill(temp_chi_sq);
            if (temp_chi_sq < best_chi_squared)
            {
                for (int i = 0; i < epsilon->GetNrows(); i++)
                {
                    double x = (x_max - x_min - bin_length) * i / (double)((n_cells_x * 2 + 1) * enlargement - 1) + x_min + bin_length / 2.0;
                    if (isnan(minimizer->X()[i]) || isnan(minimizer->Errors()[i]))
                    {
                        cout << i << "\t" << minimizer->X()[i] << "\t" << minimizer->Errors()[i] << endl;
                    }
                    best_f_vec[i] = minimizer->X()[i];
                    epsilon_best_chi_sq[i] = minimizer->Errors()[i];
                    temp_hist->Fill(x, minimizer->X()[i]);
                    temp_hist->SetBinError(i + 1, minimizer->Errors()[i]);
                }
                stringstream temp_name;
                temp_name << "debiased f with chi^2=" << to_string(temp_chi_sq);
                temp_hist->SetTitle(temp_name.str().c_str());
                temp_hist->Write(temp_name.str().c_str());
                best_chi_squared = temp_chi_sq;
                cout << "best_chi_squared=" << best_chi_squared << endl;
            }
        }
    }
    double temp_scaling = energy / (double)temp_hist->Integral();
    double mse_shower_mc_truth = 0;
    double mse_ncdfs = 0;
    TH1D *debiased_best_f_integral = new TH1D("debiased_best_f", "debiased f NCDF;x/mm;y/mm;", temp_hist->GetNbinsX(), x_min, x_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    compute_NCDF(temp_hist, debiased_best_f_integral);
    for (uint i = 0; i < temp_hist->GetNbinsX(); i++)
    {
        temp_hist->SetBinContent(i + 1, temp_hist->GetBinContent(i + 1) * temp_scaling);
        temp_hist->SetBinError(i + 1, temp_hist->GetBinErrorUp(i + 1) * temp_scaling);
        mse_shower_mc_truth += pow(hist_fine_1D->GetBinContent(i + 1) - temp_hist->GetBinContent(i + 1), 2);
        mse_ncdfs += pow(cdf_fine_1D->GetBinContent(i + 1) - debiased_best_f_integral->GetBinContent(i + 1), 2);
    }
    mse_shower_mc_truth /= temp_hist->GetNbinsX();
    chi_sq_hist->Write();
    temp_hist->Write("best fit");
    status_hist->Write();
    debiased_best_f_integral->Write("best chi sq NCDF");
    cout << "chi squared = " << best_chi_squared << endl;
    cout << "mse between shower and mc truth=" << mse_shower_mc_truth << endl;
    cout << "mse between ncdfs=" << mse_ncdfs << endl;
    hist_fine_1D->Write("MC truth");
    cdf_fine_1D->Write("MC truth NCDF");
    my_file->Close();
}

void CalorimeterShower::compute_NCDF(TH1D *in, TH1D *out)
{
    if (in->GetNbinsX() != out->GetNbinsX())
    {
        cout << "Bins do not match" << endl;
        abort();
    }
    for (uint i = 1; i <= in->GetNbinsX(); i++)
    {
        double temp = 0;
        for (uint k = 1; k <= i; k++)
        {
            temp += in->GetBinContent(k);
        }
        out->SetBinContent(i, temp / in->Integral());
    }
}

void CalorimeterShower::multiresolution(uint n_cells_x, uint n_cells_y)
{
    if (!use_toolbox && !mc_sim)
    {
        cout << "Define a toolbox file or start a (simplified) MC simulation." << endl;
        abort();
    }
    set_x_range(-(n_cells_x + 0.5) * module_dimension, (n_cells_x + 0.5) * module_dimension);
    if (n_cells_y == 0 || dimensions == 1)
    {
        multiresolution(n_cells_x);
        return 0;
    }
    set_y_range(-(n_cells_y + 0.5) * module_dimension, (n_cells_y + 0.5) * module_dimension);
    amp_results.reserve(amp_iterations);
    amp_mse = new TGraph(amp_iterations);
    amp_mse->SetTitle(";Iteration;MSE");
    TH1D *chi_sq_hist = new TH1D("chi sq", "chi squared distribution", 500, 0, 0);
    TH2D *approximation = new TH2D("approximation", "high res approximation", 100, 0, 0, 100, 0, 0);
    if (use_toolbox)
    {
        TFile *my_file = new TFile(toolbox_file_name.c_str(), "read");
        stringstream toolbox_plot;
        if (toolbox_directory_name != "")
        {
            toolbox_plot << toolbox_directory_name << "/";
        }
        toolbox_plot << toolbox_th2d_name;
        (*approximation) = (*(TH2D *)my_file->Get(toolbox_plot.str().c_str()));
        high_res_approximation = new TVectorD((2 * n_cells_x + 1) * (2 * n_cells_y + 1) * pow(enlargement, dimensions));
        epsilon = new TVectorD(high_res_approximation->GetNrows());
        int max_x, max_y, max_z;
        approximation->GetBinXYZ(approximation->GetMaximumBin(), max_x, max_y, max_z);
        uint x_lower_bound = max_x + approximation_x_shift - enlargement / 2 - n_cells_x * enlargement;
        uint x_upper_bound = max_x + approximation_x_shift + enlargement / 2 + n_cells_x * enlargement;
        uint y_lower_bound = max_y + approximation_y_shift - enlargement / 2 - n_cells_y * enlargement;
        uint y_upper_bound = max_y + approximation_y_shift + enlargement / 2 + n_cells_y * enlargement;
        for (uint bin_x = x_lower_bound; bin_x <= x_upper_bound; bin_x++)
        {
            for (uint bin_y = y_lower_bound; bin_y <= y_upper_bound; bin_y++)
            {
                uint g_index = bin_y - y_lower_bound + (bin_x - x_lower_bound) * (2 * n_cells_y + 1) * enlargement;
                (*high_res_approximation)[g_index] = approximation->GetBinContent(bin_x, bin_y);
                (*epsilon)[g_index] = approximation->GetBinErrorUp(bin_x, bin_y);
            }
        }
        my_file->Close();
    }
    stringstream output_file_name;
    output_file_name << "result_x" << enlargement << ".root";
    TFile *my_file = new TFile(output_file_name.str().c_str(), "recreate");
    TH2D *approximation_cut_out = new TH2D("approximation cut out", "Approximation (g) cut out;x/mm;y/mm", (2 * n_cells_x + 1) * enlargement, x_min, x_max, (2 * n_cells_y + 1) * enlargement, y_min, y_max);
    for (uint i = 0; i < high_res_approximation->GetNrows(); i++)
    {
        vector<unsigned int> bin_numbers(convert_index_to_bins(i, 2 * n_cells_x + 1, 2 * n_cells_y + 1));
        approximation_cut_out->SetBinContent(bin_numbers[0], bin_numbers[1], (*high_res_approximation)[i]);
        approximation_cut_out->SetBinError(bin_numbers[0], bin_numbers[1], (*epsilon)[i]);
    }
    approximation_cut_out->Write("cut out");
    cout << "Is the cut out centered?(y/n)" << endl;
    string check;
    cin >> check;
    if (check == "n" || check == "N" || check == "no" || check == "No")
    {
        abort();
    }
    TVectorD amp_solution(amp(high_res_approximation, 0, amp_iterations));
    TH2D *amp_histo = new TH2D("amp", "Solution of AMP", (n_cells_x * 2 + 1) * enlargement, x_min, x_max, (n_cells_y * 2 + 1) * enlargement, y_min, y_max);
    amp_histo->SetLineColor(kRed + 2);
    for (int i = 0; i < amp_solution.GetNrows(); i++)
    {
        unsigned int x_bin = i / ((n_cells_x * 2 + 1) * enlargement);
        unsigned int y_bin = i % ((n_cells_y * 2 + 1) * enlargement);
        amp_histo->SetBinContent(x_bin + 1, y_bin + 1, amp_solution[i]);
    }
    for (uint i = 0; i < amp_results.size(); i++)
    {
        amp_results[i].Write();
    }
    amp_results.clear();
    amp_histo->Write("amp");
    amp_mse->Write("mse");
    // debiasing
    function<double(const double *)> chisquare_data = chisquare(L, high_res_approximation, epsilon, &out_pairs_1);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(0);
    minimizer->SetTolerance(debias_tolerance);
    // minimizer->SetPrintLevel(2);
    unsigned int numbOfArguments = epsilon->GetNrows();
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments); // function of type double
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));

    double best_chi_squared = DBL_MAX;
    TH2D *temp_hist = new TH2D("best fit histogram", "debiased f fit;x/mm;y/mm;E/GeV", (n_cells_x * 2 + 1) * enlargement, x_min, x_max, (n_cells_y * 2 + 1) * enlargement, y_min, y_max);
    TVectorD best_f_vec((n_cells_x * 2 + 1) * enlargement * (n_cells_y * 2 + 1) * enlargement);
    TVectorD epsilon_best_chi_sq((n_cells_x * 2 + 1) * enlargement * (n_cells_y * 2 + 1) * enlargement);
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    if (get_variance() == 0)
    {
        cout << "0 variance. Abort." << endl;
        abort();
    }
    for (int i = 0; i < epsilon->GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    TH1D *status_hist = new TH1D("status_hist", "Minimizer status", 5, 0, 0);
    double bin_length = module_dimension / ((double)enlargement);
    bool found_fit = false;
    while (!found_fit)
    {
        while (counter < fit_attempts || (best_chi_squared > max_debias_chi_sq && found_fit))
        {
            cout << endl;
            cout << "fit attempt=" << counter + 1 << endl;
            for (int i = 0; i < epsilon->GetNrows(); i++)
            {
                stringstream name;
                name << "f" << to_string(i);
                double start_value = distributions[i](generator);
                if (start_value < 0)
                {
                    start_value = 0;
                }
                minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
            }
            minimizer->Minimize();

            if (minimizer->Status() == 4)
            {
                double temp_edm = minimizer->Edm() + 1;
                while (temp_edm > minimizer->Edm() && minimizer->Status() == 4)
                {
                    cout << "Status is 4, Edm is " << minimizer->Edm() << ". Continuing minimizing this fit attempt." << endl;
                    temp_edm = minimizer->Edm();
                    minimizer->Minimize();
                }
            }

            status_hist->Fill(minimizer->Status());
            counter++;
            if (minimizer->Status() <= 1)
            {
                vector<double> args_vec(epsilon->GetNrows(), 0);
                for (uint i = 0; i < args_vec.size(); i++)
                {
                    args_vec[i] = minimizer->X()[i];
                }
                double *args = args_vec.data();
                double temp_chi_sq = chisquare_result(args);
                cout << "chi sq = " << temp_chi_sq << endl;
                chi_sq_hist->Fill(temp_chi_sq);
                if (temp_chi_sq < best_chi_squared)
                {
                    for (int i = 0; i < epsilon->GetNrows(); i++)
                    {
                        unsigned int x_bin = i % ((n_cells_x * 2 + 1) * enlargement);
                        unsigned int y_bin = floor(i / ((n_cells_y * 2 + 1) * enlargement));
                        double x = (x_max - x_min - bin_length) * x_bin / (double)((n_cells_x * 2 + 1) * enlargement - 1) + x_min + bin_length / 2.0;
                        double y = (y_max - y_min - bin_length) * y_bin / (double)((n_cells_y * 2 + 1) * enlargement - 1) + y_min + bin_length / 2.0;
                        if (isnan(minimizer->X()[i]) || isnan(minimizer->Errors()[i]))
                        {
                            cout << i << "\t" << minimizer->X()[i] << "\t" << minimizer->Errors()[i] << endl;
                        }
                        best_f_vec[i] = minimizer->X()[i];
                        epsilon_best_chi_sq[i] = minimizer->Errors()[i];
                        temp_hist->Fill(x, y, minimizer->X()[i]);
                        temp_hist->SetBinError(x_bin + 1, y_bin + 1, minimizer->Errors()[i]);
                    }
                    stringstream temp_name;
                    temp_name << "debiased f with chi^2=" << to_string(temp_chi_sq);
                    temp_hist->SetTitle(temp_name.str().c_str());
                    temp_hist->Write(temp_name.str().c_str());
                    best_chi_squared = temp_chi_sq;
                    cout << "best_chi_squared=" << best_chi_squared << endl;
                    found_fit = true;
                }
            }
            else
            {
                cout << "Finished with status " << minimizer->Status() << "." << endl;
            }
        }
        if (!found_fit)
        {
            debias_tolerance *= debias_tolerance_increase;
            minimizer->SetTolerance(debias_tolerance);
            cout << "No fit found. Increasing tolerance by a factor of " << debias_tolerance_increase << " to " << debias_tolerance << "." << endl;
            counter = 0;
        }
    }
    TH1D *relative_error_dist = new TH1D("relative errors", "Relative Errors", 200, 0, 1);
    double temp_scaling = energy / (double)temp_hist->Integral();
    for (uint i = 0; i < temp_hist->GetNbinsX(); i++)
    {
        for (uint j = 0; j < temp_hist->GetNbinsY(); j++)
        {
            temp_hist->SetBinContent(i + 1, j + 1, temp_hist->GetBinContent(i + 1, j + 1) * temp_scaling);
            temp_hist->SetBinError(i + 1, j + 1, temp_hist->GetBinErrorUp(i + 1, j + 1) * temp_scaling);
            relative_error_dist->Fill(temp_hist->GetBinErrorUp(i + 1, j + 1) / abs(temp_hist->GetBinContent(i + 1, j + 1)));
        }
    }
    chi_sq_hist->Write();
    temp_hist->Write("best fit");
    relative_error_dist->Write();
    TH2D *debiased_best_f_integral = new TH2D("debiased_best_f", "debiased f NCDF;x/mm;y/mm;", temp_hist->GetNbinsX(), x_min, x_max, temp_hist->GetNbinsY(), y_min, y_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    compute_2D_NCDF(temp_hist, debiased_best_f_integral);
    status_hist->Write();
    debiased_best_f_integral->Write("best chi sq NCDF");
    cout << "chi squared = " << best_chi_squared << endl;
    TH1D *high_res_ncdf_x_projection = new TH1D("best_f_ncdf_x_proj", "x projection of NCDF;x/mm;y/mm", debiased_best_f_integral->GetNbinsX(), x_min, x_max);
    create_x_projection(debiased_best_f_integral, high_res_ncdf_x_projection);
    high_res_ncdf_x_projection->Write();
    construct_real_data_plots(my_file);
    my_file->Close();
}

void CalorimeterShower::calibrate_penalty()
{
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(0);
    minimizer->SetTolerance(0.1);
    // minimizer->SetPrintLevel(2);
    TF1 lambda_fnct(
        "lambda", [](double *args, double *pars)
        {
            double z = args[0];
            double result = 1 - (1 + z * z) * erf(-z) + z * exp(-z * z / 2.0) / sqrt(2.0 * TMath::Pi());
            result /= (z * z + result);
            return -result; },
        0, 1);
    ROOT::Math::Functor f = ROOT::Math::Functor(lambda_fnct, 1);
    minimizer->SetFunction(f);
    double min_function_value = DBL_MAX;
    srand(time(NULL));
    for (uint fit_attempt = 0; fit_attempt < fit_attempts_lambda; fit_attempt++)
    {
        double start_z = rand() / RAND_MAX;
        minimizer->SetLowerLimitedVariable(0, "z", start_z, 1e-5, 0.0);
        minimizer->Minimize();
        double temp_function_value = minimizer->MinValue();
        if (temp_function_value < min_function_value)
        {
            lambda = minimizer->X()[0];
            min_function_value = temp_function_value;
            cout << "update lambda to " << lambda << endl;
            cout << "min function value=" << min_function_value << endl;
        }
    }
}

unsigned int CalorimeterShower::convert_bins_to_index(unsigned int bin_x, unsigned int bin_y, unsigned int n_y)
{
    unsigned int result = (bin_y - 1) + (bin_x - 1) * n_y;
    return result;
}

vector<unsigned int> CalorimeterShower::convert_index_to_bins(unsigned int index, unsigned int n_bins_x, unsigned int n_bins_y)
{
    unsigned int y_bin = index % (n_bins_x * enlargement);
    unsigned int x_bin = floor(index / (n_bins_y * enlargement));
    vector<unsigned int> result = {x_bin + 1, y_bin + 1};
    return result;
}

void CalorimeterShower::compute_wavelet_solution()
{
    TVectorD start_value((*L_dual) * (*high_res_approximation));
    TVectorD algo3_solution(algorithm3(start_value, (*high_res_approximation), 0, 100));
    TH1D *algo3_histo = new TH1D("algorithm3", "Solution of algorithm 3;x/mm", algo3_solution.GetNrows(), x_min, x_max);
    algo3_histo->SetLineColor(kRed);
    TH1D *algo3Integral = new TH1D("algorithm3Integral", "Integral of algorithm 3;x/mm", algo3_solution.GetNrows(), x_min, x_max);
    algo3Integral->SetLineColor(kRed);
    TH1D *high_res_approximation_int = new TH1D("approx_integral", "Integral of g;x/mm", high_res_approximation->GetNrows(), x_min, x_max);
    double algo3_sum = 0;
    double high_res_approx_sum = 0;
    for (int i = 0; i < algo3_solution.GetNrows(); i++)
    {
        algo3_histo->SetBinContent(i + 1, algo3_solution[i]);
        algo3_sum += algo3_solution[i];
        algo3Integral->SetBinContent(i + 1, algo3_sum);
        high_res_approx_sum += (*high_res_approximation)[i];
        high_res_approximation_int->SetBinContent(i + 1, high_res_approx_sum);
    }
    if (write)
    {
        algo3_histo->Write();
    }

    for (int i = 0; i < algo3_solution.GetNrows(); i++)
    {
        algo3Integral->SetBinContent(i + 1, algo3Integral->GetBinContent(i + 1) / algo3_sum);
        high_res_approximation_int->SetBinContent(i + 1, high_res_approximation_int->GetBinContent(i + 1) / high_res_approx_sum);
    }

    if (write)
    {
        algo3Integral->Write();
        high_res_approximation_int->Write();
    }
    if (!use_toolbox)
    {
        for (int i = 0; i < algo3Integral->GetNbinsX(); i++)
        {
            tot_error += pow(algo3_histo->GetBinContent(i) - hist_fine_1D->GetBinContent(i), 2);
            tot_error_ncdf += pow(algo3Integral->GetBinContent(i + 1) - cdf_fine_1D->GetBinContent(i + 1), 2);
        }
        tot_error_ncdf /= algo3Integral->GetNbinsX();
        tot_error /= algo3Integral->GetNbinsX();
    }
}

void CalorimeterShower::set_real_data_showerparameters_for_shashlik(vector<double> in1, vector<double> in2)
{
    if (in1.size() == in2.size() - 1)
    {
        real_data_a_shashlik = in1;
        real_data_b_shashlik = in2;
    }
    else if (in2.size() == in1.size() - 1)
    {
        real_data_b_shashlik = in1;
        real_data_a_shashlik = in2;
    }
    else
    {
        cout << "The vector for a needs one less entry than the vector for b." << endl;
        abort();
    }
}

void CalorimeterShower::construct_real_data_plots(TFile *root_file)
{
    TF2 *ncdf(construct_ncdf_function(&real_data_a[detector_type], &real_data_b[detector_type], module_dimension, "real data ncdf"));
    ncdf->SetNpx(n_modules_x * enlargement);
    ncdf->SetNpy(n_modules_y * enlargement);
    TH2D *ncdf_histo = (TH2D *)ncdf->CreateHistogram();
    real_data_ncdf_x_projection = new TH1D("real data ncdf x projection", "X projection of real data NCDF;x/mm", ncdf_histo->GetNbinsX(), x_min, x_max);
    create_x_projection(ncdf_histo, real_data_ncdf_x_projection);
    if (root_file->IsOpen())
    {
        ncdf->Write();
        real_data_ncdf_x_projection->Write();
    }
}

string CalorimeterShower::construct_atan_sum(int n_a, vector<string> *par_names)
{
    stringstream result;
    result << "(";
    stringstream last_a;
    last_a << "(1-";
    string x_string = "(x+[" + to_string(2 * n_a + 1) + "])";
    string y_string = "(y+[" + to_string(2 * n_a + 1) + "])";
    for (uint i = 0; i < n_a; i++)
    {
        stringstream a_var;
        a_var << "[" << to_string(2 * i) << "]";
        stringstream b_var;
        b_var << "[" << to_string(2 * i + 1) << "]";
        result << a_var.str() << "*(atan(" << x_string << "/" << b_var.str() << ")+atan(" << y_string << "/" << b_var.str() << ")+atan(" << x_string << "*" << y_string << "/" << b_var.str() << "/sqrt(" << b_var.str() << "*" << b_var.str() << "+" << x_string << "*" << x_string << "+" << y_string << "*" << y_string << "))) + ";
        last_a << a_var.str();
        if (i == n_a - 1)
        {
            last_a << ")";
        }
        else
        {
            last_a << " - ";
        }
        stringstream a_name;
        a_name << "a" << to_string(i + 1);
        par_names->push_back(a_name.str());
        stringstream b_name;
        b_name << "b" << to_string(i + 1);
        par_names->push_back(b_name.str());
    }
    stringstream b_var;
    b_var << "[" << to_string(2 * n_a) << "]";
    stringstream b_name;
    b_name << "b" << to_string(1 + n_a);
    par_names->push_back(b_name.str());
    par_names->push_back("cell_width");
    result << last_a.str() << "*(atan(" << x_string << "/" << b_var.str() << ")+atan(" << y_string << "/" << b_var.str() << ")+atan(" << x_string << "*" << y_string << "/" << b_var.str() << "/sqrt(" << b_var.str() << "*" << b_var.str() << "+" << x_string << "*" << x_string << "+" << y_string << "*" << y_string << "))))/(2*TMath::Pi())+0.25";
    return result.str();
}

TF2 *CalorimeterShower::construct_ncdf_function(vector<double> *in1, vector<double> *in2, double width, string tf_name)
{
    vector<double> a, b;
    if (in1->size() == in2->size() - 1)
    {
        a = (*in1);
        b = (*in2);
    }
    else if (in2->size() == in1->size() - 1)
    {
        a = (*in2);
        b = (*in1);
    }
    else
    {
        cout << "construct_ncdf_function: Shower parameters do not fit convention." << endl;
        abort();
    }
    vector<string> par_names;
    string formula = construct_atan_sum(a.size(), &par_names);
    TF2 *ncdf = new TF2(tf_name.c_str(), formula.c_str(), x_min, x_max, y_min, y_max);
    for (uint i = 0; i < par_names.size(); i++)
    {
        ncdf->SetParName(i, par_names[i].c_str());
    }
    for (uint i = 0; i < a.size(); i++)
    {
        ncdf->SetParameter(2 * i, a[i]);
        ncdf->SetParameter(2 * i + 1, b[i]);
    }
    ncdf->SetParameter(a.size() + b.size() - 1, b.back());
    ncdf->SetParameter(a.size() + b.size(), width / 2.0);
    return ncdf;
}

void CalorimeterShower::sort_lednev_parameters(vector<double> *in1, vector<double> *in1_errors, vector<double> *in2, vector<double> *in2_errors, TMatrixD *in_cov_matrix, TH2D *out_histo)
{
    vector<double> *lednev_a;
    vector<double> *lednev_a_errors;
    vector<double> *lednev_b;
    vector<double> *lednev_b_errors;
    if (in1->size() == in2->size() - 1)
    {
        lednev_a = in1;
        lednev_a_errors = in1_errors;
        lednev_b = in2;
        lednev_b_errors = in2_errors;
    }
    else if (in2->size() == in1->size() - 1)
    {
        lednev_a = in2;
        lednev_a_errors = in2_errors;
        lednev_b = in1;
        lednev_b_errors = in1_errors;
    }
    else
    {
        cout << "sort_lednev_parameters: shower parameters do not follow convention." << endl;
        abort();
    }
    TMatrixD new_cov_matrix(numbOfArguments + 1, numbOfArguments + 1);
    for (uint i = 0; i < in_cov_matrix->GetNrows(); i++)
    {
        double temp_cov = 0;
        for (uint j = 0; j < in_cov_matrix->GetNcols(); j++)
        {
            new_cov_matrix[i][j] = (*in_cov_matrix)[i][j];
            if (j % 2 == 0 && j != in_cov_matrix->GetNcols() - 1)
            {
                temp_cov -= (*in_cov_matrix)[i][j];
            }
        }
        new_cov_matrix[i][new_cov_matrix.GetNcols() - 1] = temp_cov;
        new_cov_matrix[new_cov_matrix.GetNrows() - 1][i] = temp_cov;
    }
    double last_cov = 0;
    for (uint i = 0; i < in_cov_matrix->GetNrows() - 1; i++)
    {
        for (uint j = 0; j < in_cov_matrix->GetNcols() - 1; j++)
        {
            if (i % 2 == 0 && j % 2 == 0)
            {
                last_cov += (*in_cov_matrix)[i][j];
            }
        }
    }
    new_cov_matrix[new_cov_matrix.GetNrows() - 1][new_cov_matrix.GetNcols() - 1] = last_cov;
    for (uint i = 0; i < new_cov_matrix.GetNrows(); i++)
    {
        swap(new_cov_matrix[new_cov_matrix.GetNrows() - 1][i], new_cov_matrix[new_cov_matrix.GetNrows() - 2][i]);
    }
    for (uint i = 0; i < new_cov_matrix.GetNrows(); i++)
    {
        swap(new_cov_matrix[i][new_cov_matrix.GetNcols() - 1], new_cov_matrix[i][new_cov_matrix.GetNcols() - 2]);
    }
    double last_a = 1;
    double last_a_error = 0;
    for (uint i = 0; i < lednev_a->size(); i++)
    {
        last_a -= lednev_a->at(i);
        last_a_error += pow(lednev_a_errors->at(i), 2);
    }
    last_a_error = sqrt(last_a_error);
    lednev_a->push_back(last_a);
    lednev_a_errors->push_back(last_a_error);
    bool sorted = false;
    while (!sorted)
    {
        bool switched_once = false;
        for (uint i = 0; i < lednev_a->size(); i++)
        {
            for (uint j = i + 1; j < lednev_a->size(); j++)
            {
                if (lednev_a->at(i) < lednev_a->at(j))
                {
                    swap(lednev_a->at(i), lednev_a->at(j));
                    swap(lednev_a_errors->at(i), lednev_a_errors->at(j));
                    swap(lednev_b->at(i), lednev_b->at(j));
                    swap(lednev_b_errors->at(i), lednev_b_errors->at(j));
                    for (uint k = 0; k < new_cov_matrix.GetNrows(); k++)
                    {
                        swap(new_cov_matrix[k][i], new_cov_matrix[k][j]);
                    }
                    for (uint k = 0; k < new_cov_matrix.GetNrows(); k++)
                    {
                        swap(new_cov_matrix[i][k], new_cov_matrix[j][k]);
                    }
                    switched_once = true;
                }
            }
        }
        if (!switched_once)
        {
            sorted = true;
        }
    }
    save_correlations(lednev_a_errors, lednev_b_errors, &new_cov_matrix, out_histo);
    lednev_a->pop_back();
    lednev_a_errors->pop_back();
}

void CalorimeterShower::save_correlations(vector<double> *lednev_a_errors, vector<double> *lednev_b_errors, TMatrixD *in_cov_matrix, TH2D *out_histo)
{
    for (uint i = 0; i < in_cov_matrix->GetNrows(); i++)
    {
        for (uint j = 0; j < in_cov_matrix->GetNcols(); j++)
        {
            if (i % 2 == 0 && j % 2 == 1)
            {
                out_histo->SetBinContent(i / 2, j / 2, (*in_cov_matrix)[i][j] / sqrt(lednev_a_errors->at(i / 2) * lednev_b_errors->at(j / 2)));
            }
        }
    }
}

string CalorimeterShower::construct_shower_sum(uint n_a, vector<string> *par_names)
{
    stringstream result;
    result << "(";
    stringstream last_a;
    last_a << "(1-";
    for (uint i = 0; i < n_a; i++)
    {
        stringstream a_var;
        a_var << "[" << to_string(2 * i) << "]";
        stringstream b_var;
        b_var << "[" << to_string(2 * i + 1) << "]";
        result << a_var.str() << "*" << b_var.str() << "/pow(x*x+y*y+" << b_var.str() << "*" << b_var.str() << ",1.5)+";
        last_a << a_var.str();
        if (i == n_a - 1)
        {
            last_a << ")";
        }
        else
        {
            last_a << " - ";
        }
        stringstream a_name;
        a_name << "a" << to_string(2 * i);
        par_names->push_back(a_name.str());
        stringstream b_name;
        b_name << "b" << to_string(2 * i + 1);
        par_names->push_back(b_name.str());
    }
    stringstream b_var;
    b_var << "[" << to_string(2 * n_a) << "]";
    stringstream b_name;
    b_name << "b" << to_string(2 * n_a);
    par_names->push_back(b_name.str());
    result << last_a.str() << "*" << b_var.str() << "/pow(x*x+y*y+" << b_var.str() << "*" << b_var.str() << ",1.5))/(2*TMath::Pi())";
    return result.str();
}

TF2 *CalorimeterShower::construct_shower_function(vector<double> *in1, vector<double> *in2, string in_name)
{
    vector<double> *a, *b;
    if (in1->size() == in2->size() - 1)
    {
        a = in1;
        b = in2;
    }
    else if (in2->size() == in1->size() - 1)
    {
        a = in2;
        b = in1;
    }
    else
    {
        cout << "construct_shower_function: Shower parameters do not fit convention." << endl;
        abort();
    }
    vector<string> par_names;
    string formula = construct_shower_sum(a->size(), &par_names);
    TF2 *shower = new TF2(in_name.c_str(), formula.c_str(), x_min, x_max, y_min, y_max);
    for (uint i = 0; i < par_names.size(); i++)
    {
        shower->SetParName(i, par_names[i].c_str());
    }
    for (uint i = 0; i < a->size(); i++)
    {
        shower->SetParameter(2 * i, a->at(i));
        shower->SetParameter(2 * i + 1, b->at(i));
    }
    shower->SetParameter(2 * a->size(), b->back());
    return shower;
}

void CalorimeterShower::create_x_projection(TH2D *in_hist, TH1D *out_hist)
{
    for (uint i = 0; i < in_hist->GetNbinsX(); i++)
    {
        double temp_sum = 0;
        for (uint j = 0; j < in_hist->GetNbinsY(); j++)
        {
            temp_sum += in_hist->GetBinContent(i + 1, j + 1);
        }
        out_hist->SetBinContent(i + 1, temp_sum);
    }
    out_hist->Scale(1.0 / out_hist->Integral());
}

void CalorimeterShower::convert_TF2_NCDF_to_TH2D(TF2 *in, TH2D *out)
{
    if (in->GetNpx() != out->GetNbinsX() || out->GetNbinsY() != in->GetNpy())
    {
        TH2D *new_out = new TH2D("", out->GetTitle(), in->GetNpx(), out->GetXaxis()->GetXmin(), out->GetXaxis()->GetXmax(), in->GetNpy(), out->GetYaxis()->GetXmin(), out->GetYaxis()->GetXmax());
        (*out) = (*new_out);
    }
    double binWidthX = in->GetXaxis()->GetBinWidth(1);
    double binWidthY = in->GetYaxis()->GetBinWidth(1);
    for (uint i = 0; i < in->GetNpx(); i++)
    {
        double x = in->GetXaxis()->GetBinCenter(i + 1) + binWidthX / 2.0;
        for (uint j = 0; j < in->GetNpy(); j++)
        {
            double y = in->GetYaxis()->GetBinCenter(j + 1) + binWidthY / 2.0;
            out->SetBinContent(i + 1, j + 1, (*in)(x, y, 0, 0));
        }
    }
}

void CalorimeterShower::transform_lednev_solution(TH2D *in, TH2D *out)
{
    TVectorD solution(in->GetNbinsX() * in->GetNbinsY());
    TVectorD imageVector(in->GetNbinsX() * in->GetNbinsY());
    for (uint i = 0; i < in->GetNbinsX(); i++)
    {
        for (uint j = 0; j < in->GetNbinsY(); j++)
        {
            imageVector[convert_bins_to_index(i + 1, j + 1, in->GetNbinsY())] = in->GetBinContent(i + 1, j + 1);
        }
    }
    for (uint i = 0; i < imageVector.GetNrows(); i++)
    {
        for (vector<uint>::iterator it = out_pairs_1[i].begin(); it != out_pairs_1[i].end(); it++)
        {
            solution[i] += (*L)[i][*it] * imageVector[*it];
        }
        vector<uint> bins(convert_index_to_bins(i, in->GetNbinsX(), in->GetNbinsY()));
        out->SetBinContent(bins[0], bins[1], solution[i]);
    }
}

void CalorimeterShower::plot_red_chi_sq_contributions_per_bin(TH2D *inData, TH2D *inModel, TH2D *out)
{
    for (uint i = 0; i < inData->GetNbinsX(); i++)
    {
        for (uint j = 0; j < inData->GetNbinsY(); j++)
        {
            if (inData->GetBinErrorUp(i + 1, j + 1) == 0)
            {
                continue;
            }
            out->SetBinContent(i + 1, j + 1, pow((inData->GetBinContent(i + 1, j + 1) - inModel->GetBinContent(i + 1, j + 1)) / inData->GetBinErrorUp(i + 1, j + 1), 2));
        }
    }
}

vector<vector<uint>> CalorimeterShower::filter_indices(uint radius)
{
    vector<vector<uint>> result;
    vector<uint> radii;
    uint counter = 1;
    for (uint i = radius + 1; i < round(n_modules_x * enlargement / 2.0); i++)
    {
        for (uint j = 0; j < counter; j++)
        {
            i++;
        }
        counter++;
        radii.push_back(i);
    }
    for (uint i = 1; i <= n_modules_x * enlargement; i++)
    {
        uint r_x = abs(i - round(n_modules_x * enlargement / 2.0));
        for (uint j = 1; j <= n_modules_y * enlargement; j++)
        {
            uint r_y = abs(j - round(n_modules_y * enlargement / 2.0));
            uint r_sup = max(r_x, r_y);
            vector<uint> temp = {i, j};
            if (r_sup <= radius || find(radii.begin(), radii.end(), r_sup) != radii.end())
            {
                result.push_back(temp);
            }
        }
    }
    result.shrink_to_fit();
    return result;
}