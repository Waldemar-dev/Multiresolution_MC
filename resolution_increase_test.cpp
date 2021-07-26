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
    unsigned int mag=2;
    unsigned int n_events = 3 * pow(10, 4);
    unsigned int n_modules=7;
    unsigned int dim=2;
    vector<double> best_values={0,0,DBL_MAX};
    CalorimeterShower shower(mag,n_modules, &a, &b, Haar, id);
    unsigned int counter=0;
    if (id == Munich1 || id == Munich2 || id == Munich3){
        TGraph2D *gr=new TGraph2D();
        gr->SetTitle("squared difference;#Gamma;#sigma");
        for (double i=gamma_min;i<gamma_max;i*=10){
            for (double j=sigma_min;j<sigma_max;j*=10){
                vector<double> par={i,j};
                shower.set_par(par);
                shower.generate_shower(n_events);
                double tot_error=shower.get_tot_error();
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
    shower.write_file();
    //shower.eliminate_unphysical_entries();
    shower.generate_shower(n_events,dim);
    best_values[2]=shower.get_tot_error();
    
    cout<<endl;
    cout<<best_values[0]<<"\t"<<best_values[1]<<"\t"<<best_values[2]<<endl;
    cout<<counter<<" data points computed"<<endl;
}