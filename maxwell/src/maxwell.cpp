#include <maxwell.h>
#include <operators.h>
#include <maxwellparameters.h>
#include <iostream>
#include <cmath>

//Constructor
Maxwell::Maxwell (Domain &d, Solver &s) : ODE(2,0){
   
    if(d.getGhostPoints() < 2){
        std:: cerr << "Warning: domain has fewer ghost points than expected.";
    }
    domain = &d;
    solver = &s;
    params = nullptr;

    //Set some default parameters.
    //params = new Parameters();

    reallocateData();
}

// Destructor
Maxwell::~Maxwell(){
    //delete params;
}

//Right hand solver
void Maxwell::rhs(const Grid &grid, double **u, double **dtu){


    unsigned int nb = domain->getGhostPoints();
    double dx = grid.getSpacing();
    int nx = grid.getSize();

    // Left boundary condition, needs to be integrated in time
    dtu[U_EY][nb] =  (u[U_EY][nb+1] - u[U_EY][nb]) / (dx);
    dtu[U_BZ][nb] =  (u[U_BZ][nb+1] - u[U_BZ][nb]) / (dx);

    // Center of the grid
    for (unsigned int i = nb+1; i < nx-1-nb ; i++) {
        dtu[U_EY][i] = - (u[U_BZ][i+1] - u[U_BZ][i-1]) / (2.0*dx);
        dtu[U_BZ][i] = - (u[U_EY][i+1] - u[U_EY][i-1]) / (2.0*dx);
    }

    // Right boundary condition, needs to be integrated in time
    dtu[U_EY][nx-1-nb] = - (u[U_EY][nx-1-nb] - u[U_EY][nx-2-nb]) / (dx);
    dtu[U_BZ][nx-1-nb] = - (u[U_BZ][nx-1-nb] - u[U_BZ][nx-2-nb]) / (dx);

}

void Maxwell::setParameters(MaxwellParameters *p){
    params = p;
}




void Maxwell::initData(){
    //the center of our Gaussian

    //use our parameters
    double x0 = 0.5;
    double sigma = params->getGaussianSigma();

    //Loop through every grid and start assigning points
    for(auto it = data.begin(); it != data.end(); ++it){
        const double *x = it->getGrid().getPoints();
        unsigned int nx = it->getGrid().getSize();
        double **u = it->getData();
        for(unsigned int i = 0; i < nx; i++){
            double val = std::exp(-(x[i] - x0)*(x[i] - x0)/(sigma*sigma));
            u[U_EY][i] = val;
            u[U_BZ][i] = -val*std::sin(x[i]);
        }
    }
}
