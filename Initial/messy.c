
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double r0 = 1.0;//1.3941e16;
   double t0 = 1.0;//86400*45.;
   double rho_csm = 1.0;//1.3177e-17;

   double rc = 0.1*r0;
   double rt = 0.31654*rc;

   double day = t0/45.;
   double t_ej = 5.*day;
   double rho_ej  = 3.0307*rho_csm;

   double v0 = r0/t0;
   double Poverrho = 1e-5*v0*v0;

   double v   = 0.0;
   double X   = 0.0;
   double rho = rho_csm;

   if( r < rc ){
      rho = rho_ej*pow(rc/r,10.);
      v = r/t_ej;
      X = 1.0;
      if( r < rt ) rho = rho_ej*pow(rc/rt,10.)*rt/r;
   }
   if( r>2.*rc ){
      rho *= .0001;
   }
 
   prim[RHO] = rho;
   prim[PPP] = rho*Poverrho;
   prim[VRR] = v;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}
