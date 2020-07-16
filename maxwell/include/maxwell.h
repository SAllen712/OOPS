#ifndef MAXWELL_H
#define MAXWELL_H

#include <ode.h>
#include <maxwellparameters.h>

class Maxwell : public ODE {

    private:
    //Variable Labels
        static const unsigned int U_EY = 0;
        static const unsigned int U_BZ = 1;

        MaxwellParameters *params;

    protected:
        //virtual void applyBoundaries(bool intermediate);
        

    public:
        Maxwell(Domain& d, Solver& a);
        virtual ~Maxwell();
	virtual void rhs(const Grid& grid, double **u, double **dudt);
        inline MaxwellParameters* getParameters(){
		return params;
	};	
        void setParameters (MaxwellParameters *p);
        virtual void initData();
};

#endif
