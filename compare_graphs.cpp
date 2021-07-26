#include "CalorimeterShower.cc"
#include "ThresholdOperator.cc"
#include "DiscreteWaveletTransformation.cc"

double case1_function(double *x, double *par){
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
    if (p == 0)
    {
        result = (pow(a * a * a - 27 * c, 1 / 3.0) - a) / 3.0;
    }
    return result;
}

bool case1(double *x, double *par){
    bool result=false;
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
        result = true;
    }
    return result;
}

bool case2(double *x, double *par){
    bool result=false;
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
        if (delta <= 0 && p!=0)
        {
            result=true;
        }
        return result;
}

double case2_function(double *x, double *par){
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
    double Gamma = -q / 2.0 * sqrt(27.0 / abs(pow(p, 3)));
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
        return result;
}

bool case3(double *x, double *par){
    bool result=false;
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
    if (delta > 0 && p < 0)
        {
            result=true;
        }
        return result;
}

double case3_function(double *x, double *par){
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
    double Gamma = -q / 2.0 * sqrt(27.0 / abs(pow(p, 3)));
    if (delta > 0 && p < 0)
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
        return result;
}

bool case4(double *x, double *par){
    bool result=false;
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
    if (delta > 0 && p > 0)
        {
            result=true;
        }
        return result;
}

double case4_function(double *x, double *par){
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
    double Gamma = -q / 2.0 * sqrt(27.0 / abs(pow(p, 3)));
            double eta = (asinh(Gamma)) / 3.0;
            double thresholded_value = (2 * sqrt(3 * b - a * a) * sinh(eta) - a) / 3.0;
            result = thresholded_value;
        return result;
}

void compare_graphs(){
    double x_min=0;
    double x_max=1000;
    double Gamma=pow(10,2);
    double sigma=pow(10,0);
    double thresholdID=(double)Munich1;
    if (thresholdID==(double)Munich1){
        Gamma=pow(10,-4);
        sigma=pow(10,-3);
    } else if(thresholdID==(double) Munich2){
        Gamma=pow(10,-6);
        sigma=pow(10,-2);
    } else if (thresholdID==(double)Munich3){
        Gamma=pow(10,-5);
        sigma=pow(10,-2);
    }
    double par[3]={Gamma,sigma,thresholdID};
    TCanvas *c1=new TCanvas("c", "c", 1200, 1200);
    double divide_space=pow(10,-5);
    c1->Divide(3,2,divide_space,divide_space);
    TF1 *ls_gr=new TF1("least_squares","x",x_min,x_max);
    ls_gr->SetLineColor(kGreen);
    ls_gr->SetLineStyle(10);
    //ls_gr->Draw();
    TF1 *threshold_gr=new TF1("threshold",threshold_function,x_min,x_max,3);
    threshold_gr->SetParameters(Gamma,sigma,thresholdID);
    threshold_gr->SetLineStyle(10);
    //threshold_gr->Draw("same");   
    TGraph *gr1=new TGraph();
    gr1->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    TGraph *gr2=new TGraph();
    gr2->SetLineColor(kOrange);
    gr2->SetMarkerColor(kOrange);
    TGraph *gr3=new TGraph();
    gr3->SetLineColor(kBlue);
    gr3->SetMarkerColor(kBlue);
    TGraph *gr4=new TGraph();
    gr4->SetLineColor(kViolet);
    gr4->SetMarkerColor(kViolet);
    unsigned int counter1=0, counter2=0, counter3=0, counter4=0;
    for (uint i=0;i<pow(10,6);i++){
        double x[1]={(x_max-x_min)*i/((double)pow(10,6))+x_min};
        double y;
        if (case1(x,par)){
            y=case1_function(x,par);
            gr1->SetPoint(counter1,x[0],y);
            counter1++;
        } else if (case2(x,par)){
            y=case2_function(x,par);
            gr2->SetPoint(counter2,x[0],y);
            counter2++;
        } else if (case3(x,par)){
            y=case3_function(x,par);
            gr3->SetPoint(counter3,x[0],y);
            counter3++;
        } else if (case4(x,par)){
            y=case4_function(x,par);
            gr4->SetPoint(counter4,x[0],y);
            counter4++;
        }
    }
    c1->cd(1);
    ls_gr->Draw();
    if (counter1>0){
        gr1->Draw("same");  
        c1->cd(2);
        gr1->Draw(); 
        c1->cd(1);
    }
    if (counter2>0){
        gr2->Draw("same");
        c1->cd(3);
        gr2->SetLineColor(kOrange);
        gr2->Draw();
        c1->cd(1);
    }
    if (counter3>0){
        gr3->Draw("same");
        c1->cd(4);
        gr3->Draw();
        c1->cd(1);
    }
    if (counter4>0){
        gr4->Draw("same");
        c1->cd(5);
        gr4->Draw();
        c1->cd(0);
    }
    c1->cd(6);
    threshold_gr->Draw();
    cout<<counter1<<"\t"<<counter2<<"\t"<<counter3<<"\t"<<counter4<<endl;
}