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
        // for (uint i=0;i<b->GetNrows();i++){
        //     double scalar_product = 0;
        //     for (int j = 0; j < A->GetNcols(); j++)
        //     {
        //         scalar_product += (*A)[i][j] * args[j];
        //     }
        //     if ((*epsilon)[i]!=0){
        //         chisquare += pow((scalar_product - (*b)[i]) / (*epsilon)[i], 2);
        //     }
        // }

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

        // tbb::parallel_for(0,b->GetNrows(),[&](uint i){
        //     double scalar_product = 0;
        //     for (int j = 0; j < A->GetNcols(); j++)
        //     {
        //         scalar_product += (*A)[i][j] * args[j];
        //     }
        //     if ((*epsilon)[i]!=0){
        //         chisquare += pow((scalar_product - (*b)[i]) / (*epsilon)[i], 2);
        //     }
        // });
        return chisquare;
    };
    return out;
}