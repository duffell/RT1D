
#include "../paul.h"

static double GAMMA_LAW = 0.0;
static double RHO_FLOOR = 0.0;
static double PRE_FLOOR = 0.0;

void setHydroParams( struct domain * theDomain ){
   GAMMA_LAW = theDomain->theParList.Adiabatic_Index;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

double get_vr( double * prim ){
   return( prim[VRR] );
}

void prim2cons( double * prim , double * cons , double dV ){
   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2 = vr*vr;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);

   cons[DDD] = rho*dV;
   cons[SRR] = rho*vr*dV;
   cons[TAU] = (.5*rho*v2 + rhoe)*dV;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      cons[q] = cons[DDD]*prim[q];
   }
}

void cons2prim( double * cons , double * prim , double dV ){

   double rho = cons[DDD]/dV;
   double Sr  = cons[SRR]/dV;
   double E   = cons[TAU]/dV;

   double vr = Sr/rho;
   double v2 = vr*vr;
   double rhoe = E - .5*rho*v2;
   double gam = GAMMA_LAW;
   double Pp = (gam-1.)*rhoe;

   if( rho<RHO_FLOOR ) rho=RHO_FLOOR;
   if( Pp < PRE_FLOOR*rho ) Pp = PRE_FLOOR*rho;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[VRR] = vr;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      prim[q] = cons[q]/cons[DDD];
   }

}

void getUstar( double * prim , double * Ustar , double Sk , double Ss ){

   double rho = prim[RHO];
   double vr  = prim[VRR];
   double Pp  = prim[PPP];
   double v2  = vr*vr;

   double gam = GAMMA_LAW;

   double rhoe = Pp/(gam-1.);

   double rhostar = rho*(Sk - vr)/(Sk - Ss);
   double Pstar = Pp*(Ss - vr)/(Sk - Ss);
   double Us = rhoe*(Sk - vr)/(Sk - Ss);

   Ustar[DDD] = rhostar;
   Ustar[SRR] = rhostar*( Ss );
   Ustar[TAU] = .5*rhostar*v2 + Us + rhostar*Ss*(Ss - vr) + Pstar;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      Ustar[q] = prim[q]*Ustar[DDD];
   }

}

void flux( double * prim , double * flux ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double v2  = vr*vr;
   double gam = GAMMA_LAW;
   double rhoe = Pp/(gam-1.);
 
   flux[DDD] = rho*vr;
   flux[SRR] = rho*vr*vr + Pp;
   flux[TAU] = (.5*rho*v2 + rhoe + Pp)*vr;

   int q;
   for( q=XXX ; q<NUM_Q ; ++q ){
      flux[q] = flux[DDD]*prim[q];
   }
}

void source( double * prim , double * cons , double rp , double rm , double dVdt ){
   double Pp  = prim[PPP];
   double r  = .5*(rp+rm);
   double r2 = (rp*rp+rm*rm+rp*rm)/3.;
   cons[SRR] += 2.*Pp*(r/r2)*dVdt;
}

void source_alpha( double * prim , double * cons , double * grad_prim , double r , double dVdt ){

   double A = 0.06;
   double B = 1.0;

   double gam = GAMMA_LAW;
   double lambda = r;

   double Pp = prim[PPP];
   double rho = prim[RHO];
   double P1 = grad_prim[PPP];
   double rho1 = grad_prim[RHO];
   double g2 = -P1*rho1;
   if( g2 < 0.0 ) g2 = 0.0;
   double alpha = prim[AAA];
   double cs = sqrt(gam*fabs(Pp/rho));
   double w = cs*sqrt(alpha);

   cons[AAA] += ( A*sqrt(g2) - B*alpha*rho*w/lambda )*dVdt;
   if( cons[AAA] < 0. ) cons[AAA] = 0.;

}

double get_eta( double * prim , double * grad_prim , double r ){

   double C = 1.5e-2;
   double lambda = r;
   double gam = GAMMA_LAW;

   double cs = sqrt( gam*fabs(prim[PPP]/prim[RHO]) );

   double alpha = prim[AAA];
   if( alpha < 0.0 ) alpha = 0.0;

   double w = cs*sqrt( alpha );
   
   double eta = C*w*lambda;

   return( eta );

}

void vel( double * prim1 , double * prim2 , double * Sl , double * Sr , double * Ss ){
   
   double gam = GAMMA_LAW;

   double P1   = prim1[PPP];
   double rho1 = prim1[RHO];
   double vn1  = prim1[VRR];

   double cs1 = sqrt(fabs(gam*P1/rho1));

   double P2   = prim2[PPP];
   double rho2 = prim2[RHO];
   double vn2  = prim2[VRR];

   double cs2 = sqrt(fabs(gam*P2/rho2));

   *Ss = ( P2 - P1 + rho1*vn1*(-cs1) - rho2*vn2*cs2 )/( rho1*(-cs1) - rho2*cs2 );

   *Sr =  cs1 + vn1;
   *Sl = -cs1 + vn1;

   if( *Sr <  cs2 + vn2 ) *Sr =  cs2 + vn2;
   if( *Sl > -cs2 + vn2 ) *Sl = -cs2 + vn2;
   
}

double mindt( double * prim , double w , double r , double dr ){

   double rho = prim[RHO];
   double Pp  = prim[PPP];
   double vr  = prim[VRR];
   double gam = GAMMA_LAW;

   double cs = sqrt(fabs(gam*Pp/rho));
   double eta = get_eta( prim , NULL , r );

   double maxvr = cs + fabs( vr - w );
   double dt = dr/maxvr;
   double dt_eta = dr*dr/eta;
   if( dt > dt_eta ) dt = dt_eta;

   return( dt );

}

