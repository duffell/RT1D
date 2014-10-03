
#include "paul.h"

double get_dA( double r ){
   return( 4.*M_PI*r*r );
}

double get_moment_arm( double rp , double rm ){
   double r3 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;
   double r2 = (rp*rp + rp*rm + rm*rm)/3.;
   double r = r3/r2;
   return( r );
}

double get_dV( double rp , double rm ){
   double dr  = rp-rm;
   double r2    = (rp*rp+rm*rm+rp*rm)/3.;
   return( 4.*M_PI*r2*dr );
}
