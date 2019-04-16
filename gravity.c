
#include "paul.h"

static int grav_E_mode = 0;
static double grav_G = 0.0;

void setGravityParams( struct domain * theDomain ){
   grav_E_mode = theDomain->theParList.grav_e_mode;
   grav_G = theDomain->theParList.grav_G;
}

double get_dV( double , double );

double get_g( struct cell * c ){

   double G = grav_G;
   double rp = c->riph;
   double rm = c->riph - c->dr;
   double rc = .5*(rp+rm);

   double frac = 1.0*(rp*rp*rp - rc*rc*rc)/(rp*rp*rp - rm*rm*rm);

//   double G = 1.0; 
   double M = c->miph-frac*c->dm;
   double r = rc;

   return( -G*M/r/r );

}

void grav_src( struct cell * c , double dVdt ){

   double rho = c->prim[RHO];
   double v   = c->prim[VRR];

   double f = get_g( c );

   c->cons[SRR] += rho*f*dVdt;
   if( grav_E_mode == 0 ) c->cons[TAU] += rho*v*f*dVdt;

}

void aggregate_mass( struct domain * theDomain ){

   int rank = theDomain->rank;
   int size = theDomain->size;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;

   double M = theDomain->point_mass;
   int i;
   int imin=Ng;
   if( rank==0 ) imin = 0;

   double Mtemp = M;
   for( i=imin-1 ; i>=0 ; --i ){
      Mtemp -= theCells[i+1].dm;
      theCells[i].miph = Mtemp;
   }

   for( i=imin ; i<Nr ; ++i ){
      //if( i>=imin ) M += theCells[i].dm;
      M += theCells[i].dm;
      theCells[i].miph = M;
   }

   int imax=Nr-Ng;
   if( rank==size-1 ) imax = Nr;
   double Mtot = theCells[imax-1].miph;
   double Mtot_inner = 0.0;
   double Mrecv = 0.0;

   int nrk;
   for( nrk=0 ; nrk < size ; ++nrk ){
      if( rank==nrk ){
         Mtot_inner = Mrecv;
         Mtot += Mtot_inner;
         if( nrk<size-1 ){
            MPI_Send( &Mtot , 1 , MPI_DOUBLE , nrk+1 , 666 , MPI_COMM_WORLD );
         }
      }
      if( rank==nrk+1 ){
         MPI_Status status;
         MPI_Recv( &Mrecv , 1 , MPI_DOUBLE , nrk , 666 , MPI_COMM_WORLD , &status );
      }
   }

   for( i=0 ; i<Nr ; ++i ){
      theCells[i].miph += Mtot_inner;
   }

}

void calculate_mass( struct domain * theDomain ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;

   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp,rm;
      rp = c->riph;
      rm = 0.0;
      if( i!=0 ) rm = theCells[i-1].riph;
      double dV = get_dV( rp , rm ); 
      c->dm = c->cons[DDD];//c->prim[RHO]*dV;
   }

   aggregate_mass( theDomain );

}

void calculate_pot( struct domain * theDomain ){
   
   calculate_mass( theDomain );

   int rank=theDomain->rank;
   int size=theDomain->size;

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   double rmax = theCells[Nr-1].riph;
   double M    = theCells[Nr-1].miph;
   double pot = 0.0;
   if( rank==size-1 ) pot=-grav_G*M/rmax;

   int Ng = theDomain->Ng;
   int imax = Nr-Ng;
   if( rank==size-1 ) imax = Nr; 
   int imin = Ng;
   if( rank==0 ) imin=0;
   double Ptot = 0.0;

   for( i=imax-1 ; i>=0 ; --i ){
      struct cell * c = theCells+i;
      double dr = c->dr;
      double g = get_g( c );
      c->pot = pot + .5*g*dr;
      pot += g*dr;
      if( i==imin ) Ptot = pot;
   }

   double Ptot_outer = 0.0;
   double Precv = 0.0;

   int nrk; 
   for( nrk=size-1 ; nrk >= 0 ; --nrk ){
      if( rank==nrk ){
         Ptot_outer = Precv;
         Ptot += Ptot_outer;
         if( nrk > 0 ){
            MPI_Send( &Ptot , 1 , MPI_DOUBLE , nrk-1 , 666 , MPI_COMM_WORLD );
         }
      } 
      if( rank==nrk-1 ){
         MPI_Status status;
         MPI_Recv( &Precv , 1 , MPI_DOUBLE , nrk , 666 , MPI_COMM_WORLD , &status );
      }  
   }

   for( i=0 ; i<Nr ; ++i ){
      theCells[i].pot += Ptot_outer;
   }

}

void gravity_addsrc( struct domain * theDomain , double dt ){

   calculate_mass( theDomain );

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double rp,rm;
      rp = c->riph;
      rm = 0.0;
      if( i!=0 ) rm = theCells[i-1].riph;
      double dV = get_dV( rp , rm );
      grav_src( c , dV*dt );
   }

}


