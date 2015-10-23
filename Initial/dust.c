
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho_min = 1e-7;
   double P_min   = 1e-7;

   double rho = 1.0;
   if( r > 0.5 ) rho = rho_min;

   prim[RHO] = rho;
   prim[PPP] = P_min;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
