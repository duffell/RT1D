
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){


   double R0 = 0.1;
   double P  = 0.1;
   double a  = 0.0005;

   double drho = (.5*tanh((R0-r)/a) + .5)/pow(R0,3.);

   double rho = 1.0 + drho;
   double X   = 0.0;

   if( r < R0 ){
      X = 1.0;
   }
 
   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = 0.0;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
