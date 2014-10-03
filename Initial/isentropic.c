
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){
   double R2 = r*r;
   prim[RHO] = 1.0 + 3.0*exp(-80.*R2);
   prim[PPP] = pow(prim[RHO],5./3.);
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}
