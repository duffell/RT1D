
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   //The following is consistent with an initial
   //time t = r0/vmax = .05477

   double E = 1.0;
   double M = 1.0;

   double r0 = 0.1;
   double rho0 = 1.0;
   double Pmin = 1e-5;

   double rho = rho0;
   double v   = 0.0;
   double X   = 0.0;

   double V = 4./3.*M_PI*r0*r0*r0;
   double vmax = sqrt(10./3.*E/M);

   if( r < r0 ){
      rho += M/V;
      v = vmax*r/r0;
      X = 1.0;
   }
 
   prim[RHO] = rho;
   prim[PPP] = Pmin;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
