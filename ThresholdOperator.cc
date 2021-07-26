#include "ThresholdOperator.hh"

ThresholdOperator::ThresholdOperator(OperatorID in, Wavelets wavelet)
{
    id = in;
    waveletID = wavelet;
}

TVectorD ThresholdOperator::apply(TVectorD *x)
{
    vector<double> par={0.0,0.0};
    return apply(x, &par);
}

TVectorD ThresholdOperator::apply(TVectorD *x, TVectorD *sigmas){
    TVectorD thresholds(compute_test_threshold(x,sigmas));
    TVectorD result(x->GetNrows());
    for (uint i=0;i<result.GetNrows();i++){
        result[i]=max(0.0,(*x)[i]-thresholds[i]);
    }
    return result;
}

TVectorD ThresholdOperator::apply(TVectorD *x, vector<double> *par)
{
    TVectorD result(x->GetNrows());
    DiscreteWaveletTransformation dwt(waveletID);
    vector<double> wavelet_coefficients;
    dwt.analyse(x, &wavelet_coefficients);
    vector<double> thresholded_wavelet_coefficients;
    if (id == Donoho)
    {
        double threshold = get_donoho_threshold(&wavelet_coefficients, x->GetNrows());
        if (useGivenLambda)
        {
            threshold = lambdaDonoho;
        }
        for (uint i = 0; i < wavelet_coefficients.size(); i++)
        {
            if (abs(wavelet_coefficients[i]) - threshold <= 0)
            {
                thresholded_wavelet_coefficients.push_back(0);
            }
            else
            {
                thresholded_wavelet_coefficients.push_back((abs(wavelet_coefficients[i]) - threshold) * sgn(wavelet_coefficients[i]));
            }
        }
    }
    else if (id == Munich1 || id == Munich2 || id == Munich3)
    {
        double gamma = par->at(0);
        double n = par->at(1);
        double A = 1.0;
        double C = gamma * gamma + n;
        for (uint i = 0; i < wavelet_coefficients.size(); i++)
        {
            double B = -abs(wavelet_coefficients[i]);
            double D = -abs(wavelet_coefficients[i]) * gamma * gamma;
            double a = B / A;
            double b = C / A;
            double c = D / A;
            double p = b - a * a / 3.0;
            double q = 2 * pow(a, 3) / 27.0 - a * b / 3.0 + c;
            double delta = pow(q / 2.0, 2) + pow(p / 3.0, 3);
            if (p == 0)
            {
                thresholded_wavelet_coefficients.push_back((pow(a * a * a - 27 * c, 1 / 3.0) - a) / 3.0);
            }
            else
            {
                double Gamma = -q / 2.0 * sqrt(27.0 / abs(pow(p, 3)));
                if (delta <= 0)
                {
                    double pi = atan(1) * 4;
                    double acosTerm=acos(Gamma);
                    if (isnan(acosTerm) && Gamma>=1){
                        acosTerm=0;
                    }
                    double eta0 = (acosTerm) / 3.0;
                    double eta1 = (acosTerm + 2 * pi) / 3.0;
                    double eta2 = (acosTerm + 4 * pi) / 3.0;
                    double solution1 = (2 * sqrt(a * a - 3 * b) * cos(eta0) - a) / 3.0;
                    double solution2 = (2 * sqrt(a * a - 3 * b) * cos(eta1) - a) / 3.0;
                    double solution3 = (2 * sqrt(a * a - 3 * b) * cos(eta2) - a) / 3.0;
                    if (isnan(solution1)){abort();}
                    if (id == Munich1)
                    {
                        thresholded_wavelet_coefficients.push_back(solution1);
                    }
                    else if (id == Munich2)
                    {
                        thresholded_wavelet_coefficients.push_back(solution2);
                    }
                    else if (id == Munich3)
                    {
                        thresholded_wavelet_coefficients.push_back(solution3);
                    }
                }
                else if (delta > 0 && p < 0)
                {
                    double acoshTerm=acosh(abs(Gamma));
                    if (isnan(acoshTerm) && Gamma<=1){
                        acoshTerm=0;
                    }
                    double eta = acoshTerm / 3.0;
                    double thresholded_value=(2 * sqrt(a * a - 3 * b) * sgn(q) * cosh(eta) + a) / 3.0;
                    thresholded_wavelet_coefficients.push_back(-thresholded_value);
                }
                else if (delta > 0 && p > 0)
                {
                    double eta = (asinh(Gamma)) / 3.0;
                    double thresholded_value=(2 * sqrt(3 * b - a * a) * sinh(eta) - a) / 3.0;
                    thresholded_wavelet_coefficients.push_back(thresholded_value);
                }else {abort();}
            }
            //result[i] *= (-1);
        }
    } else if (id == Positive) {
        vector<double> x_vec;
        for (uint i=0;i<x->GetNrows();i++){
            x_vec.push_back((*x)[i]);
        }
        x_vec.shrink_to_fit();
        double threshold = get_donoho_threshold(&x_vec, x->GetNrows());
        if (useGivenLambda)
        {
            threshold = lambdaDonoho;
        }
        for (uint i=0;i<x->GetNrows();i++){
            result[i]=max(0.0,(*x)[i] - threshold);
        }
        return result;
    } else {
        cout<<"Threshold ID "<< id <<" not defined."<<endl;
        abort();
    }
    thresholded_wavelet_coefficients.shrink_to_fit();
    dwt.synthesise(&thresholded_wavelet_coefficients, &result);
    return result;
}

double ThresholdOperator::median(vector<double> in)
{
    sort(in.begin(), in.end());
    double result;
    if (in.size() % 2 == 0)
    {
        result = (in[in.size() / 2 - 1] + in[in.size() / 2]) / 2.0;
    }
    else
    {
        result = in[in.size() / 2];
    }
    return result;
}

double ThresholdOperator::mad(vector<double> in)
{
    sort(in.begin(), in.end());
    vector<double> result_vec(in.size(), 0);
    double med = median(in);
    for (uint i = 0; i < in.size(); i++)
    {
        result_vec[i] = abs(in[i] - med);
    }
    return median(result_vec);
}

double ThresholdOperator::get_donoho_threshold(vector<double> *wavelet_coefficients, unsigned int n)
{
    double sigma = mad((*wavelet_coefficients));
    double lambda = sigma * sqrt(2 * log(n)) / sqrt(n);
    return lambda;
}

double ThresholdOperator::compute_donoho_threshold(TVectorD *in)
{
    DiscreteWaveletTransformation dwt(waveletID);
    vector<double> wavelet_coefficients;
    dwt.analyse(in, &wavelet_coefficients);
    vector<double> thresholded_wavelet_coefficients;
    double threshold = get_donoho_threshold(&wavelet_coefficients, in->GetNrows());
    return threshold;
}

TVectorD ThresholdOperator::compute_test_threshold(TVectorD *in, TVectorD *sigma){
    unsigned int n=in->GetNrows();
    TVectorD result(n);
    for (uint i=0;i<n;i++){
        result[i]=(*sigma)[i] * sqrt(2 * log(n)) / sqrt(n);;
    }
    return result;
}

double ThresholdOperator::compute_threshold(TVectorD *in){
    double n=in->GetNrows();
    vector<double> in_vec;
    for (uint i=0;i<in->GetNrows();i++){
        in_vec.push_back((*in)[i]);
    }
    in_vec.shrink_to_fit();
    double sigma = mad(in_vec);
    double lambda = sigma * sqrt(2 * log(n)) / sqrt(n);
    return lambda;
}