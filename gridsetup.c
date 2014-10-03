
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"

int getN0( int drank , int dsize , int dnum ){
   int N0 = (dnum*drank)/dsize;
   return(N0);
}

void setupGrid( struct domain * theDomain ){

   int Ng = NUM_G;
   theDomain->Ng = Ng;
   int Num_R = theDomain->theParList.Num_R;
   int LogZoning = theDomain->theParList.LogZoning;

   int rank = theDomain->rank;
   int size = theDomain->size;

   double Rmin = theDomain->theParList.rmin;
   double Rmax = theDomain->theParList.rmax;

   int N0r = getN0( rank   , size , Num_R );
   int N1r = getN0( rank+1 , size , Num_R );
   if( rank != 0 ) N0r -= Ng;
   if( rank != size-1 ) N1r += Ng;
   int Nr = N1r-N0r;

   theDomain->Nr = Nr;
   theDomain->theCells = (struct cell *) malloc( Nr*sizeof(struct cell));
   printf("Rank = %d, Nr = %d\n",theDomain->rank,Nr);

   int i;

   double dx = 1./(double)Num_R;
   double x0 = (double)N0r/(double)Num_R;
   double R0 = theDomain->theParList.LogRadius;
   for( i=0 ; i<Nr ; ++i ){
      double xm = x0 + ((double)i   )*dx;
      double xp = x0 + ((double)i+1.)*dx;
      double rp,rm;
      if( LogZoning == 0 ){
         rp = Rmin + xp*(Rmax-Rmin);
         rm = Rmin + xm*(Rmax-Rmin);
      }else if( LogZoning == 1 ){
         rp = Rmin*pow(Rmax/Rmin,xp);
         rm = Rmin*pow(Rmax/Rmin,xm);
      }else{
         rp = R0*pow(Rmax/R0,xp) + Rmin-R0 + (R0-Rmin)*xp;
         rm = R0*pow(Rmax/R0,xm) + Rmin-R0 + (R0-Rmin)*xm;
      }
      theDomain->theCells[i].riph = rp;
      theDomain->theCells[i].dr   = rp - rm;
   }

}


