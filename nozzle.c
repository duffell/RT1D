
#include "paul.h"

void source_nozz( double * prim , double * cons , double rp , double rm , double t , double dVdt ){
   double r  = .5*(rp+rm);
   double R = 0.05;
   double Vol = 4./3.*M_PI*pow(R,3.);
   double T = 0.1;
   double E = 1.0;

   double Q = 0.0;
   if( r<R && t<T ) Q = E/T/Vol;

//   cons[TAU] += Q*dVdt;
}


