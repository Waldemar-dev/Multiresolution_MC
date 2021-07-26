#include <math.h>
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace TMath;

function<double(const double *)> chisquare_output(function<double(const double *)> chifunc)
{
    function<double(const double *)> out = [chifunc](const double *args) -> double { return chifunc(args); };
    return out;
}

function<double(const double *)> chisquare(TMatrixD *A, TVectorD *b, TVectorD *epsilon) //A*x=b+/-epsilon
{
    function<double(const double *)> out = [A, b, epsilon](const double *args)
        -> double {
        double chisquare = 0;
        for (uint i = 0; i < b->GetNrows(); i++)
        {
            double scalar_product = 0;
            for (uint j = 0; j < A->GetNcols(); j++)
            {
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