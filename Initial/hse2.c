
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double R = 1.0;
   double G = 1.0;
   double M = 1.0;
   
   double rho_min = 1e-8;

   int n = 5;
   double rhoc = M/R/R/R/(4.*M_PI*sqrt(3.));
   double rho = rhoc/pow( 1. + r*r/R/R/3. , 2.5 );
   double Pp = 2./3.*M_PI*G*rhoc*rhoc*R*R/pow( 1. + r*r/R/R/3. , 3. );


   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
