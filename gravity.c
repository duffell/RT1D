
#include "paul.h"

void grav_src( struct cell * c , double dt ){

   double m  = c->cons[DDD];
   double Sr = c->cons[SRR];


   double rp = c->riph;
   double rm = c->riph - c->dr;
   double rc = .5*(rp+rm);

   double frac = (rp*rp*rp - rc*rc*rc)/(rp*rp*rp - rm*rm*rm);

   double G = 1.0;
   double M = c->miph-frac*c->dm;
   double r = rc;
   //double f = menc_force( x , r );
   double f = -G*M/r/r;

   //printf("F=%e\n",m*f*dt);

   c->cons[SRR] += m*f*dt;
   c->cons[TAU] += Sr*f*dt;

}

void aggregate_mass( struct domain * theDomain ){

   int rank = theDomain->rank;
   int size = theDomain->size;
   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;

   double M=0;
   int i;
   int imin=Ng;
   if( rank==0 ) imin = 0;

   for( i=0 ; i<Nr ; ++i ){
      if( i>imin ) M += theCells[i].dm;
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
      }
      MPI_Scatter( &Mtot , 1 , MPI_DOUBLE , &Mrecv , 1 , MPI_DOUBLE , nrk , MPI_COMM_WORLD );
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
      theCells[i].dm = theCells[i].cons[DDD];
   }

   aggregate_mass( theDomain );

}

void gravity_addsrc( struct domain * theDomain , double dt ){

   calculate_mass( theDomain );

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   for( i=0 ; i<Nr ; ++i ){
      grav_src( theCells+i , dt );
   }

}


