
#include "paul.h"

double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

void plm( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double PLM = theDomain->theParList.PLM;
   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      int im = i-1;
      int ip = i+1;
      if( i==0 ) im = 0;
      if( i==Nr-1 ) ip = Nr-1;
      struct cell * c  = theCells+i;
      struct cell * cL = theCells+im;
      struct cell * cR = theCells+ip;
      double drL = cL->dr;
      double drC = c->dr;
      double drR = cR->dr;
      for( q=0 ; q<NUM_Q ; ++q ){
         double pL = cL->prim[q];
         double pC = c->prim[q];
         double pR = cR->prim[q];
         double sL = pC - pL;
         sL /= .5*( drC + drL );
         double sR = pR - pC;
         sR /= .5*( drR + drC );
         double sC = pR - pL;
         sC /= .5*( drL + drR ) + drC;
         c->grad[q] = minmod( PLM*sL , sC , PLM*sR );
      }
   }
}

