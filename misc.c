#include "paul.h"
#include <string.h>

double get_dA( double );
double get_dV( double , double );
double get_moment_arm( double , double );

double mindt( double * , double , double , double );

double getmindt( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   double dt = 1e100;
   int i;
   for( i=1 ; i<Nr-1 ; ++i ){
      int im = i-1;
      struct cell * c = theCells+i;
      double dr = c->dr;
      double r = c->riph-.5*dr;
      double wm = theCells[im].wiph;
      double wp = c->wiph;
      double w = .5*(wm+wp);
      double dt_temp = mindt( c->prim , w , r , dr );
      if( dt > dt_temp ) dt = dt_temp;
   }
   dt *= theDomain->theParList.CFL; 
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

   return( dt );
}

void initial( double * , double * );
void prim2cons( double * , double * , double );
void cons2prim( double * , double * , double );
double get_vr( double * );

void set_wcell( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr-1 ; ++i ){
      struct cell * cL = theCells+i;  
      double w = 0.0;
      if( mesh_motion ){
         struct cell * cR = theCells+i+1;
         double wL = get_vr( cL->prim );
         double wR = get_vr( cR->prim );
         w = .5*(wL + wR); 
         if( i==0 && theDomain->rank==0 ) w = 0.5*wR*(cL->riph)/(cR->riph-.5*cR->dr);//*(cR->riph - .5*cR->dr)/(cL->riph);//0.0;//2./3.*wR;
         //if( i==0 && theDomain->rank==0 ) w = 0.0;
      }
      cL->wiph = w;
   }
}

void adjust_RK_cons( struct domain * theDomain , double RK ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q] = (1.-RK)*c->cons[q] + RK*c->RKcons[q];
      }
   }
}

void move_cells( struct domain * theDomain , double RK , double dt){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      c->riph += c->wiph*dt;
   }

}

void calc_dr( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=1 ; i<Nr ; ++i ){
      int im = i-1;
      double rm = theCells[im].riph;
      double rp = theCells[i ].riph;
      double dr = rp-rm;
      theCells[i].dr = dr;
   }
   if( theDomain->rank==0 ) theCells[0].dr = theCells[0].riph;

}


void calc_prim( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , c->prim , dV );
   }

}

void plm( struct domain *);
void riemann( struct cell * , struct cell * , double , double );

void radial_flux( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   plm( theDomain );
   for( i=0 ; i<Nr-1 ; ++i ){
      struct cell * cL = theCells+i;
      struct cell * cR = theCells+i+1;
      double r = cL->riph;
      double dA = get_dA(r); 
      riemann( cL , cR , r , dA*dt );
   }

}

void source( double * , double * , double , double , double );
void source_nozz( double * , double * , double , double , double , double );
void source_alpha( double * , double * , double * , double , double );
void source_grow( double * , double * , double * , double , double );
void gravity_addsrc( struct domain * , double );

void add_source( struct domain * theDomain , double dt ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double t = theDomain->t;
   double grad[NUM_Q];

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp = c->riph;
      double rm = rp-c->dr;
      double r = get_moment_arm(rp,rm);
      double dV = get_dV(rp,rm);
      source( c->prim , c->cons , rp , rm , dV*dt );
      source_nozz( c->prim , c->cons , rp , rm , t , dV*dt );
      int inside = i>0 && i<Nr-1;
      for( q=0 ; q<NUM_Q ; ++q ){
         if( inside ){
            struct cell * cp = theCells+i+1;
            struct cell * cm = theCells+i-1;
            double dR = .5*cp->dr + c->dr + .5*cm->dr;
            grad[q] = (cp->prim[q]-cm->prim[q])/dR;
/*
            double gR = 2.*(cp->prim[q]-c->prim[q])/(cp->dr+c->dr);
            double gL = 2.*(c->prim[q]-cm->prim[q])/(c->dr+cm->dr);
            double g = gL;
            if( gL*gR < 0. ) g = 0.;
            if( fabs(gR) < fabs(g) ) g = gR;
            grad[q] = g;
*/
         }else{
            grad[q] = 0.0;
         }
      }
      source_alpha( c->prim , c->cons , grad , r , dV*dt );
      if( 0 ) source_grow( c->prim , c->cons , grad , r , dV*dt );
   }   

   int gravity_flag = theDomain->theParList.grav_flag;
   if( gravity_flag ) gravity_addsrc( theDomain , dt );

}


void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , int * rL , int * rS ){ 

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   double rmax = theCells[Nr-1].riph;
   double rmin = theCells[0].riph;
   MPI_Allreduce( MPI_IN_PLACE , &rmax , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &rmin , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );
   int Nr0 = theDomain->theParList.Num_R;
   double dr0 = rmax/(double)Nr0;
   double dx0 = log(rmax/rmin)/Nr0;
   int logscale = theDomain->theParList.LogZoning;

   double Long  = 0.0; 
   double Short = 0.0; 
   int iLong  = -1;
   int iShort = -1;

   int rank = theDomain->rank;
   int size = theDomain->size;
   int Ng   = theDomain->Ng;

   int imin = 1;
   int imax = Nr-1;
   if( rank!=0 )      imin = Ng;
   if( rank!=size-1 ) imax = Nr-Ng;

   int i;
   for( i=imin ; i<imax ; ++i ){
      struct cell * c = theCells+i;
      double dy = c->dr;
      double dx = dr0;
      if( logscale ) dx = c->riph*dx0;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }

   struct { double value ; int index ; } maxminbuf;
   maxminbuf.value = Short;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *rS = maxminbuf.index;

   maxminbuf.value = Long;
   maxminbuf.index = rank;
   MPI_Allreduce( MPI_IN_PLACE , &maxminbuf , 1 , MPI_DOUBLE_INT , MPI_MAXLOC , MPI_COMM_WORLD );
   *rL = maxminbuf.index;

   *iS = iShort;
   *iL = iLong;
   *S = Short;
   *L = Long;

}


void AMR( struct domain * theDomain ){

   double L,S;
   int iL=0;
   int iS=0;
   int rL=0;
   int rS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , &rL , &rS );
   int rank = theDomain->rank;

   //if( rank == rL ) printf("Rank %d; Long  = %e #%d\n",rank,L,iL);
   //if( rank == rS ) printf("Rank %d; Short = %e #%d\n",rank,S,iS);

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * theCells = theDomain->theCells;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;

   if( S>MaxShort && rank == rS ){
      //printf("KILL! rank = %d\n",rank);

      int iSp = iS+1;
      int iSm = iS-1;
      //Possibly shift iS backwards by 1 
      double drL = theCells[iSm].dr;
      double drR = theCells[iSp].dr;
      int imin = Ng;
      if( rank==0 ) imin = 0;
      if( drL<drR && iSm>imin ){
         --iS;
         --iSm;
         --iSp;
      }
      struct cell * c  = theCells+iS;
      struct cell * cp = theCells+iSp;

      //Remove Zone at iS+1
      c->dr   += cp->dr;
      c->riph  = cp->riph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]   += cp->cons[q];
         c->RKcons[q] += cp->RKcons[q];
      }
      double rp = c->riph;
      double rm = rp - c->dr;
      double dV = get_dV( rp , rm );
      cons2prim( c->cons , c->prim , dV );
      //Shift Memory
      int blocksize = Nr-iSp-1;
      memmove( theCells+iSp , theCells+iSp+1 , blocksize*sizeof(struct cell) );
      theDomain->Nr -= 1;
      Nr = theDomain->Nr;
      theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
      theCells = theDomain->theCells;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong && rank==rL ){

      //printf("FORGE! rank = %d\n",rank);
      theDomain->Nr += 1;
      Nr = theDomain->Nr;
      theDomain->theCells = (struct cell *) realloc( theCells , Nr*sizeof(struct cell) );
      theCells = theDomain->theCells;
      int blocksize = Nr-iL-1;
      memmove( theCells+iL+1 , theCells+iL , blocksize*sizeof(struct cell) );

      struct cell * c  = theCells+iL;
      struct cell * cp = theCells+iL+1;

      double rp = c->riph;
      double rm = rp - c->dr;
      double r0 = pow( .5*(rp*rp*rp+rm*rm*rm) , 1./3. );

      c->riph  = r0;
      c->dr    = r0-rm;
      cp->dr   = rp-r0;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         c->cons[q]    *= .5;
         c->RKcons[q]  *= .5;
         cp->cons[q]   *= .5;
         cp->RKcons[q] *= .5;
      }

      double dV = get_dV( r0 , rm );
      cons2prim( c->cons , c->prim , dV );
      dV = get_dV( rp , r0 );
      cons2prim( cp->cons , cp->prim , dV );

   }

}


