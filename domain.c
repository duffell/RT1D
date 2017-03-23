
#include "paul.h"

double get_moment_arm( double , double );
double get_dV( double , double );

int num_diagnostics( void );

void setICparams( struct domain * );
void setHydroParams( struct domain * );
void setRiemannParams( struct domain * );

void setupDomain( struct domain * theDomain ){

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->point_mass = theDomain->theParList.grav_pointmass;

   theDomain->count_steps = 0;
   theDomain->final_step = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
 
}

void initial( double * , double ); 
void prim2cons( double * , double * , double );
void cons2prim( double * , double * , double );
void restart( struct domain * ); 
void calc_dp( struct domain * );
void set_wcell( struct domain * );

void setupCells( struct domain * theDomain ){

   int i;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = &(theCells[i]);
      double rp = c->riph;
      double rm = rp - c->dr;
      c->wiph = 0.0; 
      double r = get_moment_arm( rp , rm );
      double dV = get_dV( rp , rm );
      initial( c->prim , r ); 
      prim2cons( c->prim , c->cons , dV );
      cons2prim( c->cons , c->prim , dV );
   }

}

void freeDomain( struct domain * theDomain ){
   free( theDomain->theCells );
}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , MPI_COMM_WORLD );
   if( final ) theDomain->final_step = 1;

}

void report( struct domain * );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      report( theDomain );
      if( theDomain->rank==0 ) printf("t = %.3e\n",t);
   }

   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
/*
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      //snapshot( theDomain , filename );
   }
*/
}


