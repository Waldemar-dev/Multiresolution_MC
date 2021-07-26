#include <vector>

using namespace std;

class JacobiEigenvalueAlgorithm
{
public:
    JacobiEigenvalueAlgorithm() = default;
    ~JacobiEigenvalueAlgorithm() = default;

private:
    unsigned int state;
    double s, c, t, p, y, d, r;
    vector<double> eigenvalues;
    vector<bool> changed;
    vector<unsigned int> ind;
    unsigned int maxind(TMatrixD *, unsigned int);
    vector<double> rotate(TMatrixD*);
    bool is_symmetric(TMatrixD*);
    void update(unsigned int, double);
    void compute(TMatrixD*, vector<double>*, TMatrixD*);
}