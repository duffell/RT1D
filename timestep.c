#include "paul.h"

void onestep( struct domain * , double , double , int , int );

void timestep( struct domain * theDomain , double dt ){
   
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
   }

   onestep( theDomain , 0.0 ,     dt , 1 , 0 );
   onestep( theDomain , 0.5 , 0.5*dt , 0 , 1 );

//   onestep( theDomain , 0.0 ,     dt , 1 , 1 );

   theDomain->t += dt;   
   theDomain->count_steps += 1;

}
