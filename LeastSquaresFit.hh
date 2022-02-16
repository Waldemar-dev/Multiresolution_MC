#include <math.h>
#include "TMatrixD.h"
#include "TVectorD.h"
#include <tbb/parallel_for.h>

using namespace TMath;

function<double(const double *)> chisquare_output(function<double(const double *)> chifunc)
{
    function<double(const double *)> out = [chifunc](const double *args) -> double { return chifunc(args); };
    return out;
}

function<double(const double *)> chisquare(TMatrixD *A, TVectorD *b, TVectorD *epsilon, map<uint,vector<uint> > *non_zero_indices) //A*x=b+/-epsilon
{
    function<double(const double *)> out = [A, b, epsilon, non_zero_indices](const double *args)
        -> double {
        double chisquare = 0;

        for (map<uint,vector<uint> >::iterator it=(*non_zero_indices).begin(); it!=(*non_zero_indices).end(); it++){
            double scalar_product = 0;
            uint i=it->first;
            for (vector<uint>::iterator it2=it->second.begin();it2!=it->second.end();it2++)
            {
                uint j=*it2;
                scalar_product += (*A)[i][j] * args[j];
            }
            if ((*epsilon)[i]!=0){
                chisquare += pow((scalar_product - (*b)[i]) / (*epsilon)[i], 2);
            }
        }
        return chisquare;
    };
    return out;
}

double shower_term(double a, double b, double x, double y){
  return a*atan(x/b)+a*atan(y/b)+a*atan(x*y/(b*sqrt(b*b+x*x+y*y)));
}

double cdf(double x, double y, const double *args, unsigned int numbOfArguments){
  double pi = TMath::Pi();
  vector<double> a, b;
  double last_a=1;
  for (uint i=0;i<numbOfArguments-1; i++){
    if (i % 2 == 0){
      a.push_back(args[i]);
      last_a-=args[i];
    } else {
      b.push_back(args[i]);
    }
  }
  a.push_back(last_a);
  b.push_back(args[numbOfArguments-1]);
  a.shrink_to_fit();
  b.shrink_to_fit();
  double result=0;
  for (uint i=0;i<a.size();i++){
    result+=shower_term(a[i], b[i], x, y);
  }
  result/=(2.0*pi);
  result+=0.25;
  return result;
}

function<double(const double *)> chisquare(TH2D *data, double energy, unsigned int numbOfArgs)
{
    function<double(const double *)> out = [data, energy, numbOfArgs](const double *args)
        -> double {
            double chisquare=0;
            double small_module_dimension=data->GetXaxis()->GetBinWidth(1);
            for (uint i=0;i<data->GetNbinsX();i++){
              for (uint j=0;j<data->GetNbinsY();j++){
                double x=data->GetXaxis()->GetBinCenter(i+1);
                double y=data->GetYaxis()->GetBinCenter(j+1);
                double z=data->GetBinContent(i+1,j+1);
                double eZ=data->GetBinErrorUp(i+1,j+1);
                if (eZ==0 && z==0){continue;}
                else if (eZ==0 && z!=0){cout<<"z="<<z<<endl;abort();}
                double energy_deposit=cdf(x+small_module_dimension/2.0,y+small_module_dimension/2.0,args, numbOfArgs) + cdf(x-small_module_dimension/2.0,y-small_module_dimension/2.0,args, numbOfArgs) - cdf(x-small_module_dimension/2.0,y+small_module_dimension/2.0,args, numbOfArgs) - cdf(x+small_module_dimension/2.0,y-small_module_dimension/2.0,args, numbOfArgs);
                energy_deposit*=energy;
                chisquare+=pow((energy_deposit-z)/eZ,2);
              }
            }
            return chisquare;
        };
    return out;
}

function<double(const double *)> chisquare(TH2D *data, TMatrixD *L, unsigned int numbOfArgs)
{
    function<double(const double *)> out = [data, L, numbOfArgs](const double *args)
        -> double {
            double chisquare=0;
            double small_module_dimension=data->GetXaxis()->GetBinWidth(1);
            TVectorD data_vector(data->GetNbinsX()*data->GetNbinsY());
            TVectorD data_error_vector(data_vector.GetNrows());
            TMatrixD model_f(data_vector.GetNrows(),1);
            unsigned int counter=0;
            for (uint i=0;i<data->GetNbinsX();i++){
              for (uint j=0;j<data->GetNbinsY();j++){
                double x=data->GetXaxis()->GetBinCenter(i+1);
                double y=data->GetYaxis()->GetBinCenter(j+1);
                double z=data->GetBinContent(i+1,j+1);
                double eZ=data->GetBinErrorUp(i+1,j+1);
                data_vector[counter]=z;
                data_error_vector[counter]=eZ;
                double energy_deposit=cdf(x+small_module_dimension/2.0,y+small_module_dimension/2.0,args, numbOfArgs) + cdf(x-small_module_dimension/2.0,y-small_module_dimension/2.0,args, numbOfArgs) - cdf(x-small_module_dimension/2.0,y+small_module_dimension/2.0,args, numbOfArgs) - cdf(x+small_module_dimension/2.0,y-small_module_dimension/2.0,args, numbOfArgs);
                model_f[counter][0]=energy_deposit;
                counter++;
              }
            }
            TMatrixD model_g((*L)*model_f);
            double model_g_norm=0;
            for (uint i=0;i<model_g.GetNrows();i++){
              model_g_norm+=model_g[i][0];
            }
            for (uint i=0;i<model_g.GetNrows();i++){
              if (data_error_vector[i]==0){continue;}
              model_g[i][0]=model_g[i][0]*data->Integral()/model_g_norm;
              chisquare+=pow((model_g[i][0]-data_vector[i])/data_error_vector[i],2);
            }
            return chisquare;
        };
    return out;
}