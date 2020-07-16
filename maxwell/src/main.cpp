#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>
#include <maxwell.h>
#include <polynomialinterpolator.h>
#include <iostream>
#include <maxwellparser.h>
#include <parameters.h>
#include <maxwellparameters.h>

int main (int argc, char* argv[]){
    //construct our domain and a grid to fit on it.
    if (argc < 2) {
        std::cout << "Usage: ./Maxwell <parameter file>\n";
        return 0;
    }

    MaxwellParameters params;
    MaxwellParser parser;
    parser.updateParameters(argv[1], &params); 

    char *fnames[2];
    fnames[0] = new char[16];
    fnames[1] = new char[16];
    sprintf(fnames[0], "Ey");
    sprintf(fnames[1], "Bz");


    Domain domain = Domain();
    int N = params.getGridPoints();
    double bounds[2] = {0,0};
    bounds[0] = domain.getBounds()[0];
    bounds[1] = domain.getBounds()[1];

    std::cout << "Creating grid with " << N << " points and bounds [" << bounds[0] << ", " << bounds[1] << "]" << std::endl;

    domain.addGrid(bounds, N);

    RK4 rk4 = RK4();
    PolynomialInterpolator interpolator = PolynomialInterpolator(4);

    Maxwell ode = Maxwell(domain, rk4);
    ode.setInterpolator(&interpolator);
    ode.setParameters(&params);
    ode.initData();

    double ti = 0.0;
    double tf = 5.0;
    double dt = domain.getCFL()*(domain.getGrids().begin())->getSpacing();
    unsigned int M = (tf - ti)/dt;  
    ode.dump_csv("phi00000.csv", 0, 0);
    for(unsigned int i = 0; i < M; i++){
        double t = (i + 1)*dt;
        ode.evolveStep(dt);

        char buffer[12];
        sprintf(buffer, "phi%05d.csv", i+1);
        ode.dump_csv(buffer, t, 0);
    }
    return 0;
}
