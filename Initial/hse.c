
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho_min = 1e-3;
   double R       = 0.5;
   double P0      = 2./3.*M_PI*R*R+0.1;

   double rho = 1.0;
   double P = P0 - 2./3.*M_PI*r*r;

   if( r > R ){
      rho = rho_min;
      P   = P0 - 2./3.*M_PI*R*R;
   }


   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
