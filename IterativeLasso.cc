#include "IterativeLasso.hh"

AMPLasso::AMPLasso(vector<double> in_beta){
    p=in_beta.size();
    beta=in_beta;
}

double AMPLasso::M(double epsilon){
    return 2*epsilon*log(1.0/epsilon);
}

vector<double> AMPLasso::next_beta(TMatrixD *A){
    n=y->GetNrows();
    vector<double> result;
    TMatrixD A_T();
    A_T.Transpose((*A));
    for (uint i=0;i<beta.size();i++){
        double temp_result=0;
        for (uint j=0;j<A_T.GetNcols();j++){
            temp_result+=A_T[i][j]*r[j];
        }
        temp_result/=n;
        temp_result=temp_result+beta[i];
        if (abs(temp_result) < gamma_t[i])
        {
            temp_result = 0;
        }
        else if (temp_result > gamma_t[i])
        {
            temp_result -= gamma_t[i];
        }
        else
        {
            temp_result += gamma_t[i];
        }
        result.push_back(temp_result);
    }
    result.shrink_to_fit();
    return result;
}

vector<double> AMPLasso::next_r(TMatrixD *A, TVectorD *y){
    vector<double> result;
    double b=0;
    for (uint i=0;i<beta.size();i++){
        if (beta[i]!=0){
            b+=1.0;
        }
    }
    b/=(double)n;
    for (uint i=0;i<y->GetNrows();i++){
        double temp_result=(*y)[i]+b*r[i];
        for (uint j=0;j<A->GetNcols();j++){
            temp_result-=((*A)[i][j]*beta[j]);
        }
        result.push_back(temp_result);
    }
    result.shrink_to_fit();
    return result;
}

vector<double> AMPLasso::next_tau_sq(){
    vector<double> result;
    for (uint i=0;i<beta.size();i++){
        
    }
    return result;
}