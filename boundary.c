
#include "paul.h"

double get_moment_arm( double , double );
void initial( double * , double ); 

void boundary( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int rank = theDomain->rank;
   int size = theDomain->size;

   
   if( rank == size-1 ){
      struct cell * cB = theCells+Nr-1;
      double rp = cB->riph;
      double rm = rp-cB->dr;
      double r = get_moment_arm(rp,rm);
      initial( cB->prim , r ); 
   }
/*
   if( rank == 0 ){
      struct cell * cB = theCells+0;
      struct cell * cP = theCells+1;
      cB->prim[RHO] = cP->prim[RHO];
      cB->prim[PPP] = cP->prim[PPP];
   }
*/
/*
   int ABSORB_R0 = theDomain->theParList.Absorb_BC;
         if( ABSORB_R0 ){
            struct cell * c3 = &(theCells[jk][1]);
            struct cell * c4 = &(theCells[jk][0]);
            for( q=0 ; q<NUM_Q ; ++q ){
               c4->prim[q] = c3->prim[q];
            }    
         }    
*/

}
