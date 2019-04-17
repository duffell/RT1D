
#include "../paul.h"

static double gam = 1.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
}

void initial( double * prim , double r ){

   double osc_fact = 1.0;

   double R = 1.0;
   double G = 1.0;
   double M = 1.0;
   
   r *= osc_fact;

   int n = 5;
   double rhoc = M/R/R/R/(4.*M_PI*sqrt(3.));
   double rho = rhoc/pow( 1. + r*r/R/R/3. , 2.5 );
   double Pp = 2./3.*M_PI*G*rhoc*rhoc*R*R/pow( 1. + r*r/R/R/3. , 3. );

   rho *= pow( osc_fact , 3. );
   Pp  *= pow( osc_fact , 3.*gam );

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
