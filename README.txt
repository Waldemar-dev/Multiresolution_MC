resolution_increase_test.cpp is the main script and can be compiled and executed with 
root -l resolution_increase_test.cpp++

The class CalorimeterShower simulates a shower within the ECAL and needs for the construction the magnification (mag), 
number of modules in one dimension (n_modules), CORAL shower parameters a and b, the wavelet mask ID that seems to correlate well with the input signal, 
in this case the shower profile, and lastly OperatorID, which refers to different versions of the softthreshold operator:
Donoho is the default soft threshold operator 
Positive is the right/positive half of the last one
Munich1-3 are the different solutions to the soft threshold operator with a Cauchy-prior instead of Laplacian-prior

The settings can be set with:
1. set_x_range(min, max), where min sets the lower bound of the ECAL in mm and max the upper bound
2. set_n_modules() sets the number of calorimeter cells in one dimension
3. set_n_events() sets how many random numbers from the probability distribution given by the CORAL parameters shall be drawn
4. set_par() takes a vector<double> of two doubles with the first being interpreted as Gamma, the denominator in the Cauchy distribution, 
and the second as coefficient to the negative log likelihood of the Cauchy term. This is only needed for a Munich1-3 threshold operator.
5. generate_shower(n_events, dim) takes two unsigned ints, where the first will be used for point 3 and the second can be 1 for a one-dimensional 
(only x-axis) or 2 for a two-dimensional (x- and y-axis) calorimeter. Results will be written into a root file named "mc_tests_x(mag).root".
6. histo_to_txt(TH2D*, string) writes the TH2D histogram into a textfile with the file name given as second argument. This is used for Manim.
7. set_fit_attempts() takes an unsigned int for the number of fit attempts for the debiasing of the LASSO solution
8. get_variance() return the variance derived with the Approximate Message Passing algorithm (AMP)