#include "JacobiEigenvalueAlgorithm.hh"

unsigned int JacobiEigenvalueAlgorithm::maxind(TMatrixD *A, unsigned int rowIndex)
{
    unsigned int m = rowIndex + 1;
    for (uint index = rowIndex + 2; index < A->GetNrows(); index++)
    {
        if (abs((*A)[rowIndex][index]) > abs((*A)[rowIndex][m]))
        {
            m = index;
        }
    }
    return m;
}

void JacobiEigenvalueAlgorithm::update(unsigned int k_in, double t_in)
{
    double y = eigenvalues[k_in];
    eigenvalues[k_in] = y + t_in;
    if (y == eigenvalues[k_in] && changed[k_in])
    {
        changed[k_in] = false;
        state--;
    }
    else if (y != eigenvalues[k_in] && !changed[k_in])
    {
        changed[k_in] = true;
        state++;
    }
}

vector<double> JacobiEigenvalueAlgorithm::rotate(TMatrixD *A, unsigned int k, unsigned int l, unsigned int i, unsigned int j)
{
    vector<double> result = {c * (*A)[k][l] - s * (*A)[i][j], s * (*A)[k][l] + c * (*A)[i][j]};
    return result;
}

bool JacobiEigenvalueAlgorithm::is_symmetric(TMatrixD *A)
{
    bool result = true;
    if (A->GetNrows() != A->GetNcols())
    {
        return result;
    }
    for (uint i = 0; i < A->GetNrows(); i++)
    {
        for (uint j=i+1;j<A->GetNrows();j++){
            if ((*A)[i][j]!=(*A)[j][i]){
                return false;
            }
        }
    }
    return result;
}

void JacobiEigenvalueAlgorithm::compute(TMatrixD *A, vector<double> *eigenvalues_out, TMatrixD *eigenvectors_out)
{
    if (!is_symmetric(A)){
        abort();
    }
    if (eigenvalues_out->size() > 0)
    {
        eigenvalues->erase(eigenvalues_out->begin(), eigenvalues_out->end());
    }
    eigenvectors_out->Unit();
    unsigned int n = A->GetNrows();
    eigenvectors_out->Unit();
    state = A->GetNrows();
    for (uint k = 0; k < A->GetNrows(); k++)
    {
        ind.push_back(maxind(A));
        eigenvalues.push_back((*A)[k][k]);
        changed.push_back(true);
    }
    ind.shrink_to_fit();
    eigenvalues.shrink_to_fit();
    changed.shrink_to_fit();
    while (state != 0)
    {
        unsigned int m = 1;
        for (uint k = 1; k < A->GetNrows() - 1; k++)
        {
            if (abs((*A)[k][ind[k]]) > abs((*A)[m][ind[m]]))
            {
                m = k;
            }
        }
        unsigned int k = m;
        unsigned int l = ind[m];
        p = (*A)[k][l];
        y = (eigenvalues[l] - eigenvalues[k]) / 2.0;
        d = abs(y) + sqrt(p * p + y * y);
        r = sqrt(p * p + d * d);
        c = d / r;
        s = p / r;
        t = p * p / d;
        if (y < 0)
        {
            s = -s;
            t = -t;
        }
        (*A)[k][l] = 0.0;
        update(k, -t);
        update(l, t);
        for (uint i = 0; i < k - 1; i++)
        {
            rotate(A, i, k, i, l);
            rotate(A, k, i, i, l);
            rotate(A, k, i, l, i);
        }
        for (uint i = 0; i < n; i++)
        {
            (*eigenvectors_out)[i][k] = c * (*eigenvectors_out)[i][k] - s * (*eigenvectors_out)[i][l];
            (*eigenvectors_out)[i][l] = s * (*eigenvectors_out)[i][k] + c * (*eigenvectors_out)[i][l];
        }
        ind[k] = maxind(A, k);
        ind[l] = maxind(A, l);
    }
}