#include "CalorimeterShower.hh"

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
            if (isnan(acosTerm) && Gamma >= 1)
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
            if (isnan(acoshTerm) && Gamma <= 1)
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
    if (isnan(result))
    {
        cout << "result is nan" << endl;
        abort();
    }
    return result;
}

CalorimeterShower::CalorimeterShower(unsigned int magnification, unsigned int n_modules_, vector<double> *a, vector<double> *b, Wavelets name, OperatorID id)
{
    energy = 40.0;
    wavelet = name;
    thresholdID = id;
    enlargement = magnification;
    module_dimension = 38.3; //mm
    enlarged_module_dimension = module_dimension / enlargement;
    set_n_modules(n_modules_);
    write = false;
    for (uint i = 0; i < enlargement; i++)
    {
        TMatrixD low_res_image(n_modules, 1);
        for (uint j = 0; j < n_modules; j++)
        {
            low_res_image[j][0] = 0;
        }
        low_res_images.push_back(low_res_image);
    }
    low_res_images.shrink_to_fit();
    for (uint i = 0; i < enlargement; i++)
    {
        double x = ((i - (enlargement - 1) / 2.0) * module_dimension / enlargement);
        vector<double> temp2 = {x};
        targets[i] = temp2;
    }
    file_name << "mc_tests_x" << to_string(enlargement) << ".root";
    a1_coral = a->at(0);
    a2_coral = a->at(1);
    b1_coral = b->at(0);
    b2_coral = b->at(1);
    b3_coral = b->at(2);
    // displacement = 1.0;
    // if (enlargement == 2)
    // {
    //     displacement = 1.6;
    // }
    TMatrixD H(n_modules * enlargement, n_modules * enlargement);
    TMatrixD H_dual(n_modules * enlargement, n_modules * enlargement);
    unsigned int index = module_dimension / enlarged_module_dimension;
    double factor = 1.0;
    for (uint i = 0; i <= index / 2; i++)
    {
        mask[i - 1] = 1.0 / enlargement;
        //  mask[-i] = 1.0 / (split + 1.0);
    }
    if (index / 2.0 == (int)(index / 2.0))
    {
        mask[(unsigned int)(index / 2.0) - 1] = 1.0 / 2.0 / enlargement;
        // mask[-(unsigned int)(index / 2.0)] = 1.0 / 2.0 / (split + 1);
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
    // else
    // {
    //     abort();
    // }

    for (uint i = 0; i < wavelet_masks.size(); i++)
    {
        H.Zero();
        H_dual.Zero();
        for (map<int, double>::iterator it = wavelet_masks[i].begin(); it != wavelet_masks[i].end(); it++)
        {
            if (it->first < 0)
            {
                H[0][n_modules * enlargement + it->first] = it->second;
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
                H_dual[0][n_modules * enlargement + it->first] = it->second;
            }
            else
            {
                H_dual[0][it->first] = it->second;
            }
        }
        for (uint i = 1; i < n_modules * enlargement; i++)
        {
            for (uint j = 0; j < n_modules * enlargement; j++)
            {
                unsigned int temp_index = n_modules * enlargement - 1 + j;
                if (temp_index >= n_modules * enlargement)
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
    TMatrixD result(n_modules * enlargement, n_modules * enlargement);
    for (map<int, double>::iterator it = dual_mask.begin(); it != dual_mask.end(); it++)
    {
        if (it->first < 0)
        {
            result[0][n_modules * enlargement + it->first] = it->second;
        }
        else
        {
            result[0][it->first] = it->second;
        }
    }
    for (uint i = 1; i < n_modules * enlargement; i++)
    {
        for (uint j = 0; j < n_modules * enlargement; j++)
        {
            unsigned int temp_index = n_modules * enlargement - 1 + j;
            if (temp_index >= n_modules * enlargement)
            {
                temp_index = j - 1;
            }
            result[i][j] = result[i - 1][temp_index];
        }
    }
    return result;
}

TMatrixD CalorimeterShower::compute_L()
{
    TMatrixD result(n_modules * enlargement, n_modules * enlargement);
    unsigned int index = module_dimension / enlarged_module_dimension;
    for (uint i = 0; i <= index / 2; i++)
    {
        result[0][i] = 1.0 / enlargement;
        if (n_modules * enlargement - i < result.GetNrows())
        {
            result[0][n_modules * enlargement - i] = 1.0 / enlargement;
        }
    }
    if (index / 2.0 == (int)(index / 2.0))
    {
        result[0][(unsigned int)(index / 2.0)] = 1.0 / 2.0 / enlargement;
        result[0][n_modules * enlargement - (unsigned int)(index / 2.0)] = 1.0 / 2.0 / enlargement;
    }
    for (uint i = 1; i < n_modules * enlargement; i++)
    {
        for (uint j = 0; j < n_modules * enlargement; j++)
        {
            unsigned int temp_index = n_modules * enlargement - 1 + j;
            if (temp_index >= n_modules * enlargement)
            {
                temp_index = j - 1;
            }
            result[i][j] = result[i - 1][temp_index];
        }
    }
    return result;
}

void CalorimeterShower::set_n_modules(unsigned int in)
{
    n_modules = in;
    double temp_x_max = n_modules * module_dimension / 2.0;
    set_x_range(-temp_x_max, temp_x_max);
}

void CalorimeterShower::set_x_range(double in1, double in2)
{
    x_min = min(in1, in2);
    x_max = max(in1, in2);
}

void CalorimeterShower::generate_shower(unsigned int n_events_)
{
    n_events = n_events_;
    TMatrixD L(compute_L());
    TMatrixD L_dual(compute_L_dual());
    TVectorD high_res_approximation(n_modules * enlargement);
    TVectorD high_res_control(n_modules * enlargement);
    TFile *my_file = new TFile(file_name.str().c_str(), "recreate");
    TF1 *f1 = new TF1("f1", "([0]*[1]/(pow(x,2)+pow([1],2))+[3]*[2]/(pow(x,2)+pow([2],2))+(1-[0]-[3])*[4]/(pow(x,2)+pow([4],2)))/TMath::Pi()", x_min, x_max);
    f1->SetParameter(0, a1_coral); // a1
    f1->SetParameter(1, b1_coral); // b1
    f1->SetParameter(2, b2_coral); // b2
    f1->SetParameter(3, a2_coral); // a2
    f1->SetParameter(4, b3_coral); // b3
    f1->SetLineColor(kGreen + 1);

    TF1 *f1Int = new TF1("Int", "([0]*atan(x/[1]) + [3]*atan(x/[2]) + (1-[0]-[3])*atan(x/[4]))/TMath::Pi()+0.5", x_min, x_max);
    f1Int->SetParameter(0, a1_coral); // a1
    f1Int->SetParameter(1, b1_coral); // b1
    f1Int->SetParameter(2, b2_coral); // b2
    f1Int->SetParameter(3, a2_coral); // a2
    f1Int->SetParameter(4, b3_coral); // b3
    f1Int->SetLineColor(kGreen + 1);
    if (write)
    {
        f1->Write();
        f1Int->Write();
    }

    if ((thresholdID == Munich1 || thresholdID == Munich2 || thresholdID == Munich3) && write)
    {
        double Gamma = par[0];
        double sigma = par[1];
        double x_min = 0;
        double x_max = 1000;
        TF1 *threshold_gr = new TF1("threshold_function", threshold_function, x_min, x_max, 3);
        threshold_gr->SetParameters(Gamma, sigma, (double)thresholdID);
        threshold_gr->Write();
        TF1 *ls_gr = new TF1("least_squares", "x", x_min, x_max);
        ls_gr->SetLineColor(kGreen+1);
        ls_gr->Write();
    }

    deposition_histograms.erase(deposition_histograms.begin(), deposition_histograms.end());
    for (uint i = 0; i < enlargement; i++)
    {
        string name = "Shooting at position " + to_string(i) + " (CORAL);x/mm";
        TH1D hist_temp("Deposition", name.c_str(), n_modules, x_min, x_max);
        deposition_histograms.push_back(hist_temp);
    }
    deposition_histograms.shrink_to_fit();
    TH1D *controlHist = new TH1D("control shower", "shower as observed;x/mm", n_modules, x_min, x_max);
    controlHist->SetLineColor(kBlue);
    TH1D *controlCDF = new TH1D("control CDF", "NCDF as observed;x/mm", n_modules, x_min, x_max);
    controlCDF->SetLineColor(kBlue);
    TH1D *hist_fine = new TH1D("CORAL Deposition", "High resolution energy deposition (CORAL);x/mm", n_modules * enlargement, x_min, x_max);
    hist_fine->SetLineColor(kGreen+1);
    TH1D *cdf_fine = new TH1D("NCDF", "Normalized cumulative distribution function (CORAL);x/mm", n_modules * enlargement, x_min, x_max);
    cdf_fine->SetLineColor(kGreen+1);
    TH1D *high_res_hist = new TH1D("g Deposition", "Approximate high resolution energy deposition (g);x/mm", n_modules * enlargement, x_min, x_max);
    TH1D *high_res_hist_Lf = new TH1D("L*f Deposition", "Approximate high resolution energy deposition (L*f);x/mm", n_modules * enlargement, x_min, x_max);
    for (map<unsigned int, vector<double>>::iterator it = targets.begin(); it != targets.end(); it++)
    {
        for (uint i = 0; i < n_events; i++)
        {
            double r = f1->GetRandom();
            deposition_histograms[it->first].Fill(r + it->second[0] * displacement);
            controlHist->Fill(r);
            hist_fine->Fill(r);
        }
        for (uint i = 0; i < deposition_histograms[it->first].GetNbinsX(); i++)
        {
            deposition_histograms[it->first].SetBinContent(i + 1, deposition_histograms[it->first].GetBinContent(i + 1) / (double)n_events * energy);
        }
    }
    for (uint i = 0; i < controlHist->GetNbinsX(); i++)
    {
        controlHist->SetBinContent(i + 1, controlHist->GetBinContent(i + 1) / ((double)n_events * targets.size()) * energy);
    }
    for (uint i = 0; i < hist_fine->GetNbinsX(); i++)
    {
        hist_fine->SetBinContent(i + 1, hist_fine->GetBinContent(i + 1) / ((double)n_events * targets.size()) * energy);
    }
    TVectorD coral_vec(hist_fine->GetNbinsX());

    for (uint i = 0; i < coral_vec.GetNrows(); i++)
    {
        coral_vec[i] = hist_fine->GetBinContent(i + 1);
    }
    TVectorD coral_g(L * coral_vec);
    TH1D *coral_g_hist = new TH1D("L*coral", "CORAL cross check;x/mm", coral_g.GetNrows(), x_min, x_max);
    coral_g_hist->SetLineColor(kGreen+1);
    for (uint i = 0; i < coral_g.GetNrows(); i++)
    {
        coral_g_hist->SetBinContent(i + 1, coral_g[i]);
    }
    coral_g_hist->Write();

    if (write)
    {
        controlHist->Write();
    }

    for (uint i = 1; i <= controlHist->GetNbinsX(); i++)
    {
        double tempSum = 0;
        for (uint j = 1; j <= i; j++)
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
            deposition_histograms[i].Write();
        }
    }

    //
    // rearranging low resolution images into high resolution approximation
    for (uint i = 0; i < n_modules; i++)
    {
        vector<double> bin_entries;
        for (int j = enlargement - 1; j >= 0; j--)
        {
            bin_entries.push_back(deposition_histograms[j].GetBinContent(i + 1));
            low_res_images[j][i][0] = deposition_histograms[j].GetBinContent(i + 1);
        }
        bin_entries.shrink_to_fit();
        for (uint j = 0; j < bin_entries.size(); j++)
        {
            high_res_hist->SetBinContent(i * (bin_entries.size()) + j + 1, bin_entries[j]);
            high_res_approximation[i * (bin_entries.size()) + j] = bin_entries[j]/enlargement;
        }
    }
    TVectorD epsilon(L.GetNcols());
    TGraphErrors *g_graph = new TGraphErrors(epsilon.GetNrows());
    g_graph->SetTitle("Approximation (g);x/mm;E/GeV");
    double bin_length = (x_max - x_min) / epsilon.GetNrows();
    for (uint i = 0; i < epsilon.GetNrows(); i++)
    {
        epsilon[i] = sigmaE_cell(i, &high_res_approximation)/enlargement;
        double x = (x_max - x_min - bin_length) * i / (double)epsilon.GetNrows() + x_min + bin_length / 2.0;
        g_graph->SetPoint(i, x + bin_length / 2.0, high_res_approximation[i]);
        g_graph->SetPointError(i, bin_length, epsilon[i]);
        stdev.push_back(epsilon[i]);
    }
    stdev.shrink_to_fit();
    g_graph->Write();
    // building NCDF
    double norm_constant = 0;
    for (uint j = 1; j <= n_modules * enlargement; j++)
    {
        norm_constant += hist_fine->GetBinContent(j);
    }
    for (uint i = 1; i <= n_modules * enlargement; i++)
    {
        double temp = 0;
        for (uint j = 1; j <= i; j++)
        {
            temp += hist_fine->GetBinContent(j);
        }
        cdf_fine->SetBinContent(i, temp / norm_constant);
    }
    if (write)
    {
        hist_fine->Write();
        high_res_hist->Write();
        cdf_fine->Write();
    }

    TVectorD start_value(L_dual * high_res_approximation);
    TVectorD algo3_solution(algorithm3(start_value, high_res_approximation, 0, 100));

    if (eliminate_unphysical_values)
    {
        for (uint i = 0; i < algo3_solution.GetNrows(); i++)
        {
            double temp_entry = algo3_solution[i];
            if (temp_entry < 0)
            {
                algo3_solution[i] = 0;
            }
        }
    }
    if (enlargement == 4 || enlargement == 2)
    {
        TH1D *algo3_histo = new TH1D("algorithm3", "Solution of algorithm 3;x/mm", algo3_solution.GetNrows(), x_min, x_max);
        algo3_histo->SetLineColor(kRed);
        TH1D *algo3Integral = new TH1D("algorithm3Integral", "Integral of algorithm 3;x/mm", algo3_solution.GetNrows(), x_min, x_max);
        algo3Integral->SetLineColor(kRed);
        TH1D *high_res_approximation_int = new TH1D("approx_integral", "Integral of g;x/mm", high_res_approximation.GetNrows(), x_min, x_max);
        double algo3_sum = 0;
        double high_res_approx_sum = 0;
        for (uint i = 0; i < algo3_solution.GetNrows(); i++)
        {
            algo3_histo->SetBinContent(i + 1, algo3_solution[i]);
            algo3_sum += algo3_solution[i];
            algo3Integral->SetBinContent(i + 1, algo3_sum);
            high_res_approx_sum += high_res_approximation[i];
            high_res_approximation_int->SetBinContent(i + 1, high_res_approx_sum);
        }
        //algo3_histo->SetMaximum(hist_fine->GetMaximum());
        if (write)
        {
            algo3_histo->Write();
        }

        for (uint i = 0; i < algo3_solution.GetNrows(); i++)
        {
            algo3Integral->SetBinContent(i + 1, algo3Integral->GetBinContent(i + 1) / algo3_sum);
            high_res_approximation_int->SetBinContent(i + 1, high_res_approximation_int->GetBinContent(i + 1) / high_res_approx_sum);
        }

        if (write)
        {
            algo3Integral->Write();
            high_res_approximation_int->Write();
        }
        tot_error = 0;
        tot_error_approx = 0;
        double upper_error = 0;
        double lower_error = 0;
        for (uint i = 0; i < algo3Integral->GetNbinsX(); i++)
        {
            //tot_error += pow(algo3_histo->GetBinContent(i) - hist_fine->GetBinContent(i), 2);
            tot_error += pow(algo3Integral->GetBinContent(i + 1) - cdf_fine->GetBinContent(i + 1), 2);
            tot_error_approx += pow(cdf_fine->GetBinContent(i + 1) - high_res_approximation_int->GetBinContent(i + 1), 2);
            upper_error = max(hist_fine->GetBinContent(i + 1) - algo3_solution[i], upper_error);
            lower_error = min(hist_fine->GetBinContent(i + 1) - algo3_solution[i], lower_error);
        }
        tot_error /= algo3Integral->GetNbinsX();
        tot_error_approx /= algo3Integral->GetNbinsX();
        cout << "upper_error=" << upper_error << endl;
        cout << "lower_error=" << lower_error << endl;
    }

    TVectorD amp_solution(amp(&high_res_approximation, 0, 3));
    TH1D *amp_histo = new TH1D("amp", "Solution of AMP;x/mm", amp_solution.GetNrows(), x_min, x_max);
    amp_histo->SetLineColor(kRed + 2);
    TH1D *ampIntegral = new TH1D("ampIntegral", "Integral of AMP;x/mm", amp_solution.GetNrows(), x_min, x_max);
    ampIntegral->SetLineColor(kRed + 2);
    double amp_norm = 0;
    for (uint i = 0; i < amp_solution.GetNrows(); i++)
    {
        amp_norm += amp_solution[i];
        amp_histo->SetBinContent(i + 1, amp_solution[i]);
    }
    double amp_sum = 0;
    for (uint i = 0; i < amp_solution.GetNrows(); i++)
    {
        amp_sum += amp_solution[i];
        ampIntegral->SetBinContent(i + 1, amp_sum / amp_norm);
    }
    amp_histo->Write();
    ampIntegral->Write();
    cout << "variance=" << get_variance() << endl;
    //fit for f
    function<double(const double *)> chisquare_data = chisquare(&L, &high_res_approximation, &epsilon);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetTolerance(0.1);
    //minimizer->SetPrintLevel(1);
    unsigned int numbOfArguments = epsilon.GetNrows();
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments); // function of type double
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));

    double min_chi_squared = DBL_MAX;
    double best_chi_squared = DBL_MAX;
    TGraphErrors *best_f_fit = new TGraphErrors(epsilon.GetNrows());
    best_f_fit->SetTitle("debiased f fit;x/mm");
    TH1D *debiased_best_f_integral = new TH1D("debiased_best_f", "debiased f integral;x/mm", n_modules * enlargement, x_min, x_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    TGraphErrors *min_f_fit = new TGraphErrors(epsilon.GetNrows());
    min_f_fit->SetTitle("debiased f fit;x/mm");
    TH1D *debiased_min_f_integral = new TH1D("debiased_min_f", "debiased f integral;x/mm", n_modules * enlargement, x_min, x_max);
    debiased_min_f_integral->SetLineColor(kBlack);
    TVectorD epsilon_min_chi_sq(n_modules * enlargement);
    TVectorD epsilon_best_chi_sq(n_modules * enlargement);
    double mse_min_f_int=0;
    double mse_best_f_int=0;
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    for (uint i = 0; i < epsilon.GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    while (counter < fit_attempts)
    {
        cout << "counter=" << counter << endl;
        for (uint i = 0; i < epsilon.GetNrows(); i++)
        {
            stringstream name;
            name << "f" << to_string(i);
            double start_value = distributions[i](generator);
            if (start_value < 0)
            {
                start_value = 0;
            }
            minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
            //minimizer->SetVariable(i, name.str().c_str(), start_value, 1e-5);
        }
        minimizer->Minimize();
        counter++;
        if (minimizer->Status() <= 1)
        {
            double bin_length = (x_max - x_min) / epsilon.GetNrows();
            vector<double> args_vec(n_modules * enlargement, 0);
            for (uint i = 0; i < n_modules * enlargement; i++)
            {
                args_vec[i] = minimizer->X()[i];
            }
            double *args = args_vec.data();
            double temp_chi_sq = chisquare_result(args) / (n_modules * enlargement);
            if (temp_chi_sq < min_chi_squared)
            {
                double debiased_sum = 0;
                for (uint i = 0; i < epsilon.GetNrows(); i++)
                {
                    double x = (x_max - x_min - bin_length) * i / (double)epsilon.GetNrows() + x_min + bin_length / 2.0;
                    min_f_fit->SetPoint(i, x + bin_length / 2.0, minimizer->X()[i]);
                    min_f_fit->SetPointError(i, bin_length, minimizer->Errors()[i]);
                    epsilon_min_chi_sq[i]=minimizer->Errors()[i];
                    debiased_sum += minimizer->X()[i];
                    debiased_min_f_integral->SetBinContent(i + 1, debiased_sum);
                    
                }
                
                mse_min_f_int=0;
                for (uint i = 0; i < debiased_min_f_integral->GetNbinsX(); i++)
                {
                    debiased_min_f_integral->SetBinContent(i + 1, debiased_min_f_integral->GetBinContent(i + 1) / debiased_sum);
                    mse_min_f_int+=pow(cdf_fine->GetBinContent(i+1)-debiased_min_f_integral->GetBinContent(i + 1),2);
                }
                mse_min_f_int/=epsilon.GetNrows();
                min_chi_squared=temp_chi_sq;
            }
            if (abs(temp_chi_sq - 1) < abs(best_chi_squared - 1))
            {
                double debiased_sum = 0;
                for (uint i = 0; i < epsilon.GetNrows(); i++)
                {
                    double x = (x_max - x_min - bin_length) * i / (double)epsilon.GetNrows() + x_min + bin_length / 2.0;
                    best_f_fit->SetPoint(i, x + bin_length / 2.0, minimizer->X()[i]);
                    best_f_fit->SetPointError(i, bin_length, minimizer->Errors()[i]);
                    epsilon_best_chi_sq[i]=minimizer->Errors()[i];
                    debiased_sum += minimizer->X()[i];
                    debiased_best_f_integral->SetBinContent(i + 1, debiased_sum);
                }
                mse_best_f_int=0;
                for (uint i = 0; i < debiased_best_f_integral->GetNbinsX(); i++)
                {
                    debiased_best_f_integral->SetBinContent(i + 1, debiased_best_f_integral->GetBinContent(i + 1) / debiased_sum);
                    mse_best_f_int+=pow(cdf_fine->GetBinContent(i+1)-debiased_best_f_integral->GetBinContent(i + 1),2);
                }
                mse_best_f_int/=epsilon.GetNrows();
                best_chi_squared = temp_chi_sq;
                cout << "best_chi_squared=" << best_chi_squared << endl;
            }
        }
    }
    cout<<"lowest red chi^2="<<min_chi_squared<<endl;
    min_f_fit->Write("min red chi sq");
    debiased_min_f_integral->Write("min red chi sq int");
    cout << "best red chi^2=" << best_chi_squared << endl;
    best_f_fit->Write("best red chi sq");
    debiased_best_f_integral->Write("best red chi sq int");
    cout<<"mean squared error of min f int="<<mse_min_f_int<<endl;
    cout<<"mean squared error of best f int="<<mse_best_f_int<<endl;
    TGraphErrors *min_red_chi_g=new TGraphErrors();
    TGraphErrors *best_red_chi_g=new TGraphErrors();
    for (uint i=0;i<L.GetNrows();i++){
        double bin_length = (x_max - x_min) / L.GetNrows();
        double x = (x_max - x_min - bin_length) * i / (double)L.GetNrows() + x_min + bin_length / 2.0;
        double min_chi_sq_value=0;
        double best_chi_sq_value=0;
        double min_chi_sq_error=0;
        double best_chi_sq_error=0;
        for (uint j=0;j<L.GetNcols();j++){
            min_chi_sq_value+=L[i][j]*min_f_fit->GetY()[j];
            min_chi_sq_error+=pow(L[i][j]*epsilon_min_chi_sq[j],2);
            best_chi_sq_value+=L[i][j]*best_f_fit->GetY()[j];
            best_chi_sq_error+=pow(L[i][j]*epsilon_best_chi_sq[j],2);
        }
        min_chi_sq_error=sqrt(min_chi_sq_error);
        best_chi_sq_error=sqrt(best_chi_sq_error);
        min_red_chi_g->SetPoint(i,x+bin_length/2.0,min_chi_sq_value);
        min_red_chi_g->SetPointError(i, bin_length / 2.0,min_chi_sq_error);
        best_red_chi_g->SetPoint(i,x+bin_length/2.0,best_chi_sq_value);
        best_red_chi_g->SetPointError(i, bin_length / 2.0,best_chi_sq_error);
    }
    min_red_chi_g->Write("min_red_chi_sq_g");
    best_red_chi_g->Write("best_red_chi_sq_g");
}

TVectorD CalorimeterShower::T(TVectorD arg, unsigned int J)
{
    TMatrixD L(compute_L());
    TMatrixD L_dual(compute_L_dual());
    TMatrixD L_dual_J = L_dual;
    TMatrixD L_J(L);
    TVectorD sigmas(stdev.size());
    if (stdev.size() < 1)
    {
        cout << "stdev size is 0" << endl;
        abort();
    }
    for (uint i = 0; i < stdev.size(); i++)
    {
        sigmas[i] = stdev[i];
    }
    for (uint i = 1; i < J; i++)
    {
        L_dual_J *= L_dual;
        L_J *= L;
    }
    TVectorD result(L_dual_J * L_J * arg);
    ThresholdOperator thresholdOperator(thresholdID, wavelet);
    for (uint j = 0; j < J - 1; j++)
    {
        TMatrixD L_d_j(L_dual.GetNrows(), L_dual.GetNcols());
        L_d_j.UnitMatrix();
        TMatrixD L_j(L.GetNrows(), L.GetNcols());
        L_j.UnitMatrix();
        for (uint i = 0; i < j; i++)
        {
            L_d_j *= L_dual;
            L_j *= L;
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
                // double penalty=34;
                // double sigma=10000;
                // vector<double> par={penalty, sigma};
                temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &par);
            }
            else if (thresholdID == Test)
            {
                temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &sigmas);
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
    TMatrixD L(compute_L());
    TMatrixD L_dual(compute_L_dual());
    TVectorD new_f = L_dual * g;
    TVectorD sigmas(stdev.size());
    for (uint i = 0; i < stdev.size(); i++)
    {
        sigmas[i] = stdev[i];
    }
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
        TMatrixD L_dual_J = L_dual;
        TMatrixD L_J = L;
        for (uint i = 1; i < J; i++)
        {
            L_dual_J *= L_dual;
            L_J *= L;
        }

        TVectorD result(L_dual_J * L_J * new_f);
        for (uint j = 0; j < J - 1; j++)
        {
            TMatrixD L_d_j(L_dual.GetNrows(), L_dual.GetNcols());
            L_d_j.UnitMatrix();
            TMatrixD L_j(L.GetNrows(), L.GetNcols());
            L_j.UnitMatrix();
            for (uint i = 0; i < j; i++)
            {
                L_d_j *= L_dual;
                L_j *= L;
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
                    // double penalty=0.5;
                    // double sigma=1;
                    // vector<double> par={penalty, sigma};
                    temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &par);
                }
                else if (thresholdID == Test)
                {
                    temp += H_dual_vec[i] * thresholdOperator.apply(&tempArg, &sigmas);
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

double CalorimeterShower::sigmaE_cell(unsigned int i, TVectorD *in)
{
    // double norm = 0;
    // for (uint j = 0; j < in->GetNrows(); j++)
    // {
    //     norm += (*in)[j];
    // }
    double result = sqrt((*in)[i] / energy) * sigmaE();
    return result;
}

double CalorimeterShower::sigmaE()
{
    //double c1 = 0.0114570;
    //double c2 = 0.00145212;
    double c1 = 0.15;//in CORAL: EC02P1__ParamEnergyError
    double c2 = 0.015;
    double c3 = 0.05;

    double result = sqrt(c1 * c1 * energy + c2 * c2 * energy * energy + c3 * c3);
    return result;
}

TVectorD CalorimeterShower::amp(TVectorD *g, TVectorD *x, TVectorD *last_x, TVectorD *z, double gamma_threshold, unsigned int counter, unsigned int counter_max)
{
    TMatrixD L(pow(n_modules*enlargement,dimensions),pow(n_modules*enlargement,dimensions));
    if (dimensions==1){
        L=compute_L();
    }
    else if (dimensions==2){
        L=compute_L2();
    }
    TMatrixD L_T(L.GetNcols(), L.GetNrows());
    cout << "counter=" << counter << endl;
    L_T.Transpose(L);
    ThresholdOperator thresholdOperator(Positive, Haar);
    lambda = gamma * thresholdOperator.compute_threshold(x);
    double next_gamma_threshold = 0;
    double threshold = lambda + gamma_threshold;
    TVectorD next_z((*g) - L * (*x));
    for (uint b = 0; b < next_z.GetNrows(); b++)
    {
        for (uint c = 0; c < next_z.GetNrows(); c++)
        {
            for (uint j = 0; j < last_x->GetNrows(); j++)
            {
                double temp_arg = 0;
                for (uint d = 0; d < next_z.GetNrows(); d++)
                {
                    temp_arg += L[d][j] * (*z)[d] + L[d][j] * L[d][j] * (*last_x)[j];
                }
                if (temp_arg > gamma_threshold)
                {
                    next_z[b] = next_z[b] + L[b][j] * L[c][j] * (*z)[c];
                    next_gamma_threshold += L[b][j] * L[c][j];
                }
            }
        }
    }
    next_gamma_threshold *= threshold;
    cout << "next_gamma_threshold=" << next_gamma_threshold << endl;
    variance = next_gamma_threshold;
    TVectorD next_x(L_T * next_z);
    for (uint i = 0; i < next_x.GetNrows(); i++)
    {
        double temp = 0;
        for (uint a = 0; a < next_z.GetNrows(); a++)
        {
            temp += L[a][i] * next_z[a] + L[a][i] * L[a][i] * (*last_x)[i];
        }
        if (temp > variance)
        {
            next_x[i] = temp - variance;
        }
        else
        {
            next_x[i] = 0;
        }
    }
    cout<<"next_x="<<endl;
    next_x.Print();
    if (counter < counter_max)
    {
        TVectorD result(amp(g, &next_x, x, &next_z, next_gamma_threshold, counter + 1, counter_max));
        return result;
    }
    return next_x;
}

TVectorD CalorimeterShower::amp(TVectorD *g, unsigned int counter, unsigned int counter_max)
{
    TMatrixD L(compute_L());
    TVectorD x(pow(L.GetNcols(),dimensions));
    x.Zero();
    TVectorD last_x(x.GetNrows());
    last_x.Zero();
    ThresholdOperator temp_threshold(thresholdID, wavelet);
    double gamma_threshold = temp_threshold.compute_threshold(g);
    return amp(g, &x, &last_x, g, gamma_threshold, counter, counter_max);
}

void CalorimeterShower::generate_2D_shower(unsigned int n_events_){
    n_events = n_events_;
    write=true;
    TMatrixD L2(compute_L2());
    TFile *my_file = new TFile(file_name.str().c_str(), "recreate");
    TF2 *f2 = new TF2("f2", "([0]*[1]/pow(x*x+y*y+[1]*[1],1.5)+[2]*[3]/pow(x*x+y*y+[3]*[3],1.5)+(1-[0]-[2])*[4]/pow(x*x+y*y+[4]*[4],1.5))/(2*TMath::Pi())",x_min,x_max,x_min,x_max);
    f2->SetParNames("a1","b1","a2","b2","b3");
    f2->SetParameter(0,a1_coral);
    f2->SetParameter(1,b1_coral);
    f2->SetParameter(2,a2_coral);
    f2->SetParameter(3,b2_coral);
    f2->SetParameter(4,b3_coral);
    TF2 *f2Int=new TF2("f2_NCDF","([0]*(atan(x/[1])+atan(y/[1])+atan(x*y/[1]/sqrt([1]*[1]+x*x+y*y))) +[2]*(atan(x/[3])+atan(y/[3])+atan(x*y/[3]/sqrt([3]*[3]+x*x+y*y))) +(1-[0]-[2])*(atan(x/[4])+atan(y/[4])+atan(x*y/[4]/sqrt([4]*[4]+x*x+y*y))))/(2*TMath::Pi())+0.25", x_min,x_max,x_min,x_max);
    f2Int->SetParNames("a1","b1","a2","b2","b3");
    f2Int->SetParameter(0,a1_coral);
    f2Int->SetParameter(1,b1_coral);
    f2Int->SetParameter(2,a2_coral);
    f2Int->SetParameter(3,b2_coral);
    f2Int->SetParameter(4,b3_coral);

    if (write){
        f2->Write();
        f2Int->Write();
    }
    deposition_histograms_2d.erase(deposition_histograms_2d.begin(), deposition_histograms_2d.end());
    for (uint i = 0; i < enlargement*enlargement; i++)
    {
        string name = "Shooting at position " + to_string(i) + " (CORAL);x/mm";
        TH2D hist_temp("Deposition", name.c_str(), n_modules, x_min, x_max, n_modules, x_min,x_max);
        deposition_histograms_2d.push_back(hist_temp);
    }
    deposition_histograms_2d.shrink_to_fit();
    unsigned int target_counter=0;
    for (uint i = 0; i < enlargement; i++)
    {
        double x = ((i - (enlargement - 1) / 2.0) * module_dimension / enlargement);
        for (uint j=0;j<enlargement;j++){
            double y = ((j - (enlargement - 1) / 2.0) * module_dimension / enlargement);
            vector<double> temp2 = {x,y};
            targets[target_counter] = temp2;
            target_counter++;
        }
    }
    TH2D *controlHist = new TH2D("control shower", "shower as observed;x/mm", n_modules, x_min, x_max, n_modules, x_min, x_max);
    controlHist->SetLineColor(kBlue);
    TH2D *controlCDF = new TH2D("control CDF", "NCDF as observed;x/mm", n_modules, x_min, x_max, n_modules, x_min, x_max);
    controlCDF->SetLineColor(kBlue);
    TH2D *hist_fine = new TH2D("CORAL Deposition", "High resolution energy deposition (CORAL);x/mm", n_modules * enlargement, x_min, x_max, n_modules * enlargement, x_min, x_max);
    hist_fine->SetLineColor(kGreen+1);
    TH2D *cdf_fine = new TH2D("NCDF", "Normalized cumulative distribution function (CORAL);x/mm", n_modules * enlargement, x_min, x_max, n_modules * enlargement, x_min, x_max);
    cdf_fine->SetLineColor(kGreen+1);
    TH2D *high_res_hist = new TH2D("g Deposition", "Approximate high resolution energy deposition (g);x/mm", n_modules * enlargement, x_min, x_max, n_modules * enlargement, x_min, x_max);
    TH2D *high_res_hist_Lf = new TH2D("L*f Deposition", "Approximate high resolution energy deposition (L*f);x/mm", n_modules * enlargement, x_min, x_max, n_modules * enlargement, x_min, x_max);
    for (map<unsigned int, vector<double>>::iterator it = targets.begin(); it != targets.end(); it++)
    {
        for (uint i = 0; i < n_events; i++)
        {
            double x,y;
            f2->GetRandom2(x,y);
            deposition_histograms_2d[it->first].Fill(x + it->second[0], y+it->second[1]);
            controlHist->Fill(x,y,1);
            hist_fine->Fill(x,y,1);
        }
        for (uint i = 0; i < deposition_histograms_2d[it->first].GetNbinsX(); i++)
        {
            for (uint j=0;j<deposition_histograms_2d[it->first].GetNbinsY();j++){
                deposition_histograms_2d[it->first].SetBinContent(i + 1, j+1, deposition_histograms_2d[it->first].GetBinContent(i + 1, j+1) / (double)n_events * energy);
            }
        }
        stringstream file_name;
        file_name<<"deposition"<<to_string(it->first)<<".dat";
        histo_to_txt(&(deposition_histograms_2d[it->first]),file_name.str());
    }
    for (uint i = 0; i < controlHist->GetNbinsX(); i++)
    {
        for (uint j=0;j<controlHist->GetNbinsY();j++){
            controlHist->SetBinContent(i + 1,j+1, controlHist->GetBinContent(i + 1,j+1) / ((double)n_events) * energy);
        }
    }
    for (uint i = 0; i < hist_fine->GetNbinsX(); i++)
    {
        for (uint j=0;j<hist_fine->GetNbinsX();j++){
            hist_fine->SetBinContent(i + 1,j+1, hist_fine->GetBinContent(i + 1,j+1) / ((double)n_events) * energy);
        }
    }
    TVectorD coral_vec(hist_fine->GetNbinsX()*hist_fine->GetNbinsY());

    unsigned int coral_vec_counter=0;
    for (uint i = 0; i < hist_fine->GetNbinsX(); i++)
    {
        for (uint j=0;j<hist_fine->GetNbinsY();j++){
            coral_vec[coral_vec_counter] = hist_fine->GetBinContent(i + 1, j+1);
            coral_vec_counter++;
        }
        
    }
    TVectorD coral_g(L2 * coral_vec);
    TH1D *coral_g_hist = new TH1D("L*coral", "CORAL cross check;x/mm", coral_g.GetNrows(), 0, coral_g.GetNrows());
    coral_g_hist->SetLineColor(kGreen+1);
    for (uint i = 0; i < coral_g.GetNrows(); i++)
    {
        coral_g_hist->SetBinContent(i + 1, coral_g[i]);
    }
    coral_g_hist->Write();

    if (write)
    {
        controlHist->Write();
    }
    for (uint i = 1; i <= controlHist->GetNbinsX(); i++)
    {
        for (uint j=1;j<=controlHist->GetNbinsY();j++){
            double tempSum = 0;
            for (uint k=1;k<=i;k++){
                for (uint l=1;l<=j;l++){
                    tempSum += controlHist->GetBinContent(k,l);
                }
            }
            controlCDF->SetBinContent(i,j, tempSum / (double)controlHist->Integral());
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
    for (uint k=0;k<deposition_histograms_2d.size();k++){
        unsigned int x_offset=floor(k /(double) enlargement);
        unsigned int y_offset=k%enlargement;
        for (uint i = 0; i < n_modules; i++)
        {
            double x_bin=enlargement*(1+i)-x_offset;
            for (uint j=0;j<n_modules;j++){       
                double y_bin=enlargement*(j+1)-y_offset;
                high_res_hist->SetBinContent( x_bin+1,y_bin+1, deposition_histograms_2d[k].GetBinContent(i+1,j+1));
           }
        }
    }
    high_res_hist->Write();
    // rewriting high resolution approximation in vector notation
    TVectorD high_res_approximation(L2.GetNcols());
    for (uint i=0;i<high_res_hist->GetNbinsX();i++){
        for (uint j=0;j<high_res_hist->GetNbinsY();j++){
            high_res_approximation[i*high_res_hist->GetNbinsY()+j]=high_res_hist->GetBinContent(i+1,j+1);
        }
    }
    TVectorD epsilon(L2.GetNcols());
    TGraphErrors *g_graph = new TGraphErrors(epsilon.GetNrows());
    g_graph->SetTitle("Approximation (g);x;E/GeV");
    double bin_length = (x_max - x_min) / (double)(n_modules*enlargement);
    unsigned int non_zero_epsilon_entries=0;
    for (uint i = 0; i < epsilon.GetNrows(); i++)
    {
        epsilon[i] = sigmaE_cell(i, &high_res_approximation)/enlargement;
        if (epsilon[i]!=0){
            non_zero_epsilon_entries++;
        }
        //double x = (x_max - x_min - bin_length) * i / (double)(n_modules*enlargement) + x_min + bin_length / 2.0;
        g_graph->SetPoint(i, i, high_res_approximation[i]);
        g_graph->SetPointError(i, 0, epsilon[i]);
        stdev.push_back(epsilon[i]);
    }
    stdev.shrink_to_fit();
    g_graph->Write();
    for (uint i = 1; i <= n_modules * enlargement; i++)
    {
        for (uint j=1;j<=n_modules * enlargement; j++){
            double temp = 0;
            for (uint k=1;k<=i;k++){
                for (uint l=1;l<=j;l++){
                    temp += hist_fine->GetBinContent(k,l);
                }
            }
            cdf_fine->SetBinContent(i,j, temp / hist_fine->Integral());
        }
    }
    if (write)
    {
        hist_fine->Write();
        cdf_fine->Write();
    }
    //Intermediate step: Apply AMP to get LASSO solution and error
    TVectorD amp_solution(amp(&high_res_approximation, 0,1));
    TH2D *amp_histo = new TH2D("amp", "Solution of AMP", n_modules*enlargement, x_min, x_max, n_modules*enlargement, x_min, x_max);
    amp_histo->SetLineColor(kRed + 2);
    TH2D *ampIntegral = new TH2D("ampIntegral", "Integral of AMP", n_modules*enlargement, x_min, x_max, n_modules*enlargement, x_min, x_max);
    ampIntegral->SetLineColor(kRed + 2);
    for (uint i = 0; i < amp_solution.GetNrows(); i++)
    {
        unsigned int x_bin=i / (n_modules*enlargement);
        unsigned int y_bin=i%(n_modules*enlargement);
        amp_histo->SetBinContent(x_bin,y_bin, amp_solution[i]);
    }
    amp_histo->Write();
    cout << "variance=" << get_variance() << endl;
    //fit for f
    function<double(const double *)> chisquare_data = chisquare(&L2, &high_res_approximation, &epsilon);
    function<double(const double *)> chisquare_result = chisquare_output(chisquare_data);
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    minimizer->SetMaxFunctionCalls(1000000);
    minimizer->SetTolerance(0.1);
    minimizer->SetPrintLevel(2);
    unsigned int numbOfArguments = epsilon.GetNrows();
    ROOT::Math::Functor f = ROOT::Math::Functor(chisquare_data, numbOfArguments); // function of type double
    minimizer->SetFunction(f);
    unsigned int counter = 0;
    srand(time(NULL));

    double best_chi_squared = DBL_MAX;
    TGraph2DErrors *best_f_fit = new TGraph2DErrors(epsilon.GetNrows());
    best_f_fit->SetTitle("debiased f fit;x/mm;y/mm;E/GeV");
    TVectorD best_f_vec(pow(n_modules*enlargement,2));
    TH2D *debiased_best_f_integral = new TH2D("debiased_best_f", "debiased f NCDF;x/mm;y/mm;E/GeV", n_modules * enlargement, x_min, x_max, n_modules * enlargement, x_min, x_max);
    debiased_best_f_integral->SetLineColor(kBlack);
    TVectorD epsilon_best_chi_sq(pow(n_modules * enlargement,dimensions));
    double mse_best_f_int=0;
    default_random_engine generator;
    vector<normal_distribution<double>> distributions;
    for (uint i = 0; i < epsilon.GetNrows(); i++)
    {
        normal_distribution<double> distribution(amp_solution[i], sqrt(get_variance()));
        distributions.push_back(distribution);
    }
    distributions.shrink_to_fit();
    while (counter < fit_attempts)
    {
        cout << "counter=" << counter << endl;
        for (uint i = 0; i < epsilon.GetNrows(); i++)
        {
            stringstream name;
            name << "f" << to_string(i);
            double start_value = distributions[i](generator);
            if (start_value < 0)
            {
                start_value = 0;
            }
            minimizer->SetLowerLimitedVariable(i, name.str().c_str(), start_value, 1e-5, 0.0);
            //minimizer->SetVariable(i, name.str().c_str(), start_value, 1e-5);
        }
        minimizer->Minimize();
        counter++;
        if (minimizer->Status() <= 1)
        {
            double bin_length = (x_max - x_min) / (n_modules*enlargement);
            vector<double> args_vec(n_modules * enlargement, 0);
            for (uint i = 0; i < n_modules * enlargement; i++)
            {
                args_vec[i] = minimizer->X()[i];
            }
            double *args = args_vec.data();
            double temp_chi_sq = chisquare_result(args) / pow(non_zero_epsilon_entries,dimensions);
            if (abs(temp_chi_sq - 1) < abs(best_chi_squared - 1))
            {
                double debiased_sum = 0;
                mse_best_f_int=0;
                for (uint i = 0; i < epsilon.GetNrows(); i++)
                {
                    unsigned int x_bin=i%(n_modules*enlargement);
                    unsigned int y_bin=i/(n_modules*enlargement);
                    double x = (x_max - x_min - bin_length) * x_bin / (double)(n_modules*enlargement) + x_min + bin_length / 2.0;
                    double y = (x_max - x_min - bin_length) * y_bin / (double)(n_modules*enlargement) + x_min + bin_length / 2.0;
                    best_f_fit->SetPoint(i, x + bin_length / 2.0, y+bin_length/2.0, minimizer->X()[i]);
                    best_f_fit->SetPointError(i, bin_length/2.0, bin_length/2.0, minimizer->Errors()[i]);
                    best_f_vec[i]=minimizer->X()[i];
                    epsilon_best_chi_sq[i]=minimizer->Errors()[i];
                    debiased_sum+=minimizer->X()[i];
                }
                for (uint i=0;i<debiased_best_f_integral->GetNbinsX();i++){
                    for (uint j=0;j<debiased_best_f_integral->GetNbinsY();j++){
                        double temp_sum=0;
                        for (uint i2=0;i2<i;i2++){
                            for (uint j2=0;j2<j;j2++){
                                unsigned int temp_index=i2*debiased_best_f_integral->GetNbinsY()+j2;
                                temp_sum+=minimizer->X()[temp_index];
                            }
                        }
                        debiased_best_f_integral->SetBinContent(i + 1, j+1, temp_sum/debiased_sum);
                        mse_best_f_int+=pow(cdf_fine->GetBinContent(i+1)-debiased_best_f_integral->GetBinContent(i + 1, j + 1),2);
                    }
                }
                mse_best_f_int/=epsilon.GetNrows();
                best_chi_squared = temp_chi_sq;
                cout << "best_chi_squared=" << best_chi_squared << endl;
            }
        }
    }
    best_f_fit->Write("best red chi sq f");
    debiased_best_f_integral->Write("best red chi sq NCDF");
    cout<<"mean squared error of best f int="<<mse_best_f_int<<endl;
    TGraphErrors *best_red_chi_g=new TGraphErrors();
    TH1D *deviations_hist=new TH1D("dev", "distribution of nominator with 0 denominator", 100, 0,0);
    for (uint i=0;i<L2.GetNrows();i++){
        double bin_length = (x_max - x_min) / L2.GetNrows();
        double x = (x_max - x_min - bin_length) * i / (double)L2.GetNrows() + x_min + bin_length / 2.0;
        double best_chi_sq_value=0;
        double best_chi_sq_error=0;
        for (uint j=0;j<L2.GetNcols();j++){
            best_chi_sq_value+=L2[i][j]*best_f_fit->GetY()[j];
            best_chi_sq_error+=pow(L2[i][j]*epsilon_best_chi_sq[j],2);
        }
        best_chi_sq_error=sqrt(best_chi_sq_error);
        best_red_chi_g->SetPoint(i,x+bin_length/2.0,best_chi_sq_value);
        best_red_chi_g->SetPointError(i, bin_length / 2.0,best_chi_sq_error);
        if (epsilon[i]==0){
            deviations_hist->Fill(best_chi_sq_value);
        }
    }
    best_red_chi_g->Write("best red chi sq g");
    deviations_hist->Write();
}

TMatrixD CalorimeterShower::compute_L2(){
    TMatrixD Lx(compute_L());
    TMatrixD Ly(Lx);
    TMatrixD result(Lx.GetNrows()*Ly.GetNrows(),Lx.GetNcols()*Ly.GetNcols());
    for (uint i_x=0;i_x<Lx.GetNrows();i_x++){
        for (uint i_y=0;i_y<Ly.GetNrows();i_y++){
            for (uint j_x=0;j_x<Lx.GetNcols();j_x++){
                for (uint j_y=0;j_y<Ly.GetNcols();j_y++){
                    result[i_x*Ly.GetNrows()+i_y][j_x*Ly.GetNcols()+j_y]=Lx[i_x][j_x]*Ly[i_y][j_y];
                }
            }
        }
    }
    return result;
}

void CalorimeterShower::generate_shower(unsigned int n_events_, unsigned int dim){
    dimensions=dim;
    if (dim==1){
        generate_shower(n_events_);
    } else if (dim==2){
        generate_2D_shower(n_events_);
    }
}

void CalorimeterShower::histo_to_txt(TH2D *in, string file_name){
    ofstream file;
    file.open(file_name);
    for (uint i=0;i<in->GetNbinsX();i++){
        for (uint j=0;j<in->GetNbinsY();j++){
            file<<in->GetBinContent(i+1,j+1)<<"\t";
        }
        file<<endl;
    }
    file.close();
}