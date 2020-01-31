
#include "../paul.h"

static double gam = 1.0;

void setICparams( struct domain * theDomain ){
   gam = theDomain->theParList.Adiabatic_Index;
}

double get_fit_vel( double r ){

   double a = 0.330363;
   double b = 3.84388;
   double c = 0.00475402;
   double d = 5.13514;
   double e = 4.75756;

   double vel = ( a*r + c*pow(r,b) )*exp( -pow(r/d,e) );

   return(vel);

}

void initial( double * prim , double r ){

//   double osc_fact = 1.0;

   double R = 1.0;///sqrt(2.*M_PI);
   double G = 1.0;
   double M = 1.0;///sqrt(2.*M_PI);
   double rho_min = 1e-7*M/R/R/R;
   double k = 2.;

//   r *= osc_fact;

   double rhoc = M/R/R/R/4./M_PI/M_PI;
   double rho = rhoc*sin(r/R)/(r/R) + rho_min;
   if( r > M_PI*R ) rho = rho_min*pow(M_PI*R/r,k);
   double Pp = 2.*M_PI*G*rho*rho*R*R;
   double Pext = G*M*rho_min/(1.+k)/(M_PI*R);
   Pp += Pext;
//   if( Pp < Pext ) Pp = Pext;
   if( r > M_PI*R ) Pp = G*M*rho_min/(1.+k)/(M_PI*R)*pow(M_PI*R/r,1.+k);

   Pp += 1e-7*G*M*M/R/R/R/R;

   double A = 0.0;
   Pp *= 1. + A*exp( -pow(r/.05/R,2.)/2./M_PI );

//   rho *= pow( osc_fact , 3. );
//   Pp  *= pow( osc_fact , 3.*gam );

   double v = 0.01*sqrt(G*M/R)*get_fit_vel( r/R );

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = v;
   prim[XXX] = 0.0;
   prim[AAA] = 0.0;

}
