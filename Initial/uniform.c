
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){
   prim[RHO] = 1.0;
   prim[PPP] = 1.0;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;
}
