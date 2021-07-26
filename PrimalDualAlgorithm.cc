#include "PrimalDualAlgorithm.hh"

void PrimalDualAlgorithm::set_tuning(double in_t1, double in_t2, double in_t3){
    t1=in_t1;
    t2=in_t2;
    t3=in_t3;
}

void PrimalDualAlgorithm::set_start_values(TVectorD in_xi, TVectorD in_x, TVectorD in_x_bar){
    xi.push_back(in_xi);
    x.push_back(in_x);
    x_bar.push_back(in_x_bar);
}

void PrimalDualAlgorithm::set_operator(function<TVectorD(TVectorD,double)> in_func){
    threshold_operator=in_func;
}

TVectorD PrimalDualAlgorithm::compute(TVectorD *y, TMatrixD *A, TMatrixD *A_dual, double lambda){
    TVectorD result(x[0].GetNrows());
    iterations=1000;
    gamma=1/lambda;
    if (xi.size()<1 || x.size()<1 || x_bar.size()<1){
        abort();
    }
    for (uint i=0;i<iterations;i++){
        TVectorD xi_temp(gamma/(gamma+t1)*(xi[0]+t1*((*A)*x_bar[0]-y)));
        xi.push_back(xi_temp);
        TVectorD x_temp(threshold_operator(x[0]-t2*(*A_dual)*xi[1],t2));
        x.push_back(x_temp);
        TVectorD x_bar_temp(x[1]+t3*(x[1]-x[0]));
        x_bar.push_back(x_bar_temp);

        xi.erase(xi.begin());
        x.erase(x.begin());
        x_bar.erase(x_bar.begin());
    }
    for (uint i=0;i<result.GetNrows();i++){
        result[i]=x[1][i];
    }
    return result;
}