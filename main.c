
#include "paul.h"

int mpiSetup( struct domain * , int , char *[] );

void setupGrid( struct domain * );

void timestep( struct domain * , double );
void setupCells( struct domain * );
void exchangeData( struct domain * );
void set_wcell( struct domain * );
double getmindt( struct domain * );

void read_par_file( struct domain * );

void setupDomain( struct domain * );
void freeDomain( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );

void start_clock( struct domain * );
void generate_log( struct domain * );

int main( int argc , char * argv[] ){
 
   MPI_Init(&argc,&argv);
   struct domain theDomain = {0};

   start_clock( &theDomain ); 
   read_par_file( &theDomain );
 
   int error = mpiSetup(&theDomain,argc,argv);
   if( error==1 ) return(0);

   if(theDomain.rank==0) remove("abort");

   setupGrid( &theDomain );   
   setupDomain( &theDomain );

   setupCells( &theDomain );

   exchangeData( &theDomain );

   if( theDomain.rank==0 ){
      FILE * rFile = fopen("report.dat","w");
      fclose(rFile);
   }

   while( !(theDomain.final_step) ){

      set_wcell( &theDomain );
      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );

   }

   possiblyOutput( &theDomain , 1 );

   generate_log( &theDomain );
   MPI_Barrier(MPI_COMM_WORLD);
   freeDomain( &theDomain );
   MPI_Finalize();

   return(0);

}

