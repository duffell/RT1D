
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double r ){

   double rho0 = 1.0;
   double rho1 = 1.0;

   double R    = 0.5;
   double k    = 4.0;

   double PA = 2./3.*M_PI*rho0*rho0*(R*R-r*r);
   double PB = 4./3.*M_PI*rho0*rho1*R*R/(k+1.);
   double PC = 4.*M_PI*rho1*rho1*R*R/(k-3.)/(k+1.);
   double PD =-4.*M_PI*rho1*rho1*R*R/(k-3.)/(2.*k-2.);

   double rho = rho0;
   double P = PA + PB + PC + PD;

   if( r > R ){
      rho = rho1*pow(R/r,k);
      P   = PB*pow(R/r,1.+k) + PC*pow( R/r , 1.+k ) + PD*pow( R/r , 2.*k-2. ) ;
   }

//   if( r < .25 ) P += 10.;

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[VRR] = 0.0;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
