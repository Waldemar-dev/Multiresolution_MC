#include "CalorimeterShower.cc"
#include "ThresholdOperator.cc"
#include "DiscreteWaveletTransformation.cc"

void resolution_increase_test(){
    vector<double> a={1.2231,-0.1421};
    vector<double> b={6.2825,30.4875,42.4935};
    OperatorID id=Positive;
    double gamma_max=pow(10,4);
    double sigma_max=pow(10,7);
    double gamma_min=pow(10,-6);
    double sigma_min=pow(10,-6);
    unsigned int mag=5;
    unsigned int n_events = 5 * pow(10, 4);
    unsigned int n_modules=5;
    unsigned int dim=2;
    vector<double> best_values={0,0,DBL_MAX};
    CalorimeterShower shower(mag,n_modules, n_modules, dim, GAMSRH);
    // shower.limit_variables();
    shower.set_number_of_parameter_pairs(2);
    shower.add_initial_limits_for_a(1.7,2.0);
    // shower.add_initial_limits_for_a(-0.1,-0.2);
    // shower.add_initial_limits_for_a(-4,4);
    // shower.add_initial_limits_for_a(-4,4);
    shower.add_initial_limits_for_b(3.8,4.5);
    shower.add_initial_limits_for_b(25,35);
    // shower.add_initial_limits_for_b(50,120);
    // shower.add_initial_limits_for_b(50,100);
    // shower.add_initial_limits_for_b(100,200);
    //shower.fix_b2(40);
    unsigned int counter=0;
    if (id == Munich1 || id == Munich2 || id == Munich3){
        TGraph2D *gr=new TGraph2D();
        gr->SetTitle("squared difference;#Gamma;#sigma");
        for (double i=gamma_min;i<gamma_max;i*=10){
            for (double j=sigma_min;j<sigma_max;j*=10){
                vector<double> par={i,j};
                shower.set_par(par);
                shower.generate_shower(n_events);
                double tot_error=shower.get_mse();
                gr->SetPoint(counter,log10(i),log10(j),tot_error);
                if (tot_error<best_values[2]){
                    best_values[0]=i;
                    best_values[1]=j;
                    best_values[2]=tot_error;
                }
                counter++;
            }
        }
        //gr->Draw("surf1");
        vector<double> par={best_values[0],best_values[1]};
        shower.set_par(par);
    }
    shower.calibrate_penalty();
    shower.write_file();
    //shower.eliminate_unphysical_entries();
    // shower.generate_shower(n_events,dim);
    
    // shower.reduce_fit_range_x(2);
    // shower.reduce_fit_range_y(2);
    // shower.set_direct_fit(true);
    shower.shift_approximation_in_x(0);
    shower.shift_approximation_in_y(0);
    shower.read_plot("my_shower_cell1988_40GeV_x5_moved.root", "highResApprox");
    // shower.multiresolution(2,2);
    shower.lednev_fit("best fit;1");
    best_values[2]=shower.get_mse();
    // shower.write_stats_into_file("amp_stats.dat");
    
    // cout<<endl;
    // cout<<best_values[0]<<"\t"<<best_values[1]<<"\t"<<best_values[2]<<endl;
    // cout<<counter<<" data points computed"<<endl;
}