
#include "paul.h"

void AMR( struct domain * ); 

void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );
void move_cells( struct domain * , double , double );
void calc_dr( struct domain * );
void calc_prim( struct domain * );

void radial_flux( struct domain * , double );
void add_source( struct domain * , double );

void boundary( struct domain * );
void exchangeData( struct domain * );
void calculate_mass( struct domain * );

void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step ){
   
   adjust_RK_cons( theDomain , RK );

   radial_flux( theDomain , dt );
   add_source( theDomain , dt );

   if( first_step ) move_cells( theDomain , RK , dt );
   calc_dr( theDomain );

   calc_prim( theDomain );

   if( last_step ){
      AMR( theDomain );
   }
   boundary( theDomain );
   exchangeData( theDomain );

}

