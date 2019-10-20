
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


//   double rc = pow( rp*rm*( rp*rp + rp*rm + rm*rm )/3. , .25 );

   double frac = 1.0*(rp*rp*rp - rc*rc*rc)/(rp*rp*rp - rm*rm*rm);

//   double G = 1.0; 
   double M = c->miph-frac*c->dm;
   double r = rc;
   if( r==0 ) r = .5*(rp+rm);


   double r2_3 = (rp*rp + rm*rm + rp*rm)/3.;
   double r3_4 = (rp*rp*rp + rp*rp*rm + rp*rm*rm + rm*rm*rm)/4.;
   double r4_5 = ( pow(rp,4.) + pow(rp,3.)*rm + rp*rp*rm*rm + rp*pow(rm,3.) + pow(rm,4.) )/5.;

   double rhot = c->dm/(rp*rp*rp - rm*rm*rm);
   double mt = c->miph - c->dm - rm*rm*rm*rhot;

   double term1 = mt*mt/rp/rm;
   if( rm == 0 ) term1 = 0.0;
   double term2 = (rp+rm)*mt*rhot;
   double term3 = r4_5*rhot*rhot;
   
   double e = (mt + r3_4*rhot)/r2_3;
   double e2 = (term1+term2+term3)/r2_3;

 
//   return( -G*M/r/r );
   return( -G*sqrt(e2) );
//   return( -G*e );

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
      theCells[i].miph = Mtemp;
      Mtemp -= theCells[i+1].dm;
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
      //double rp,rm;
      //rp = c->riph;
      //rm = rp-c->dr;
      //if( i!=0 ) rm = theCells[i-1].riph;
      //double dV = get_dV( rp , rm ); 
      c->dm = c->cons[DDD];//c->prim[RHO]*dV;
   }

   aggregate_mass( theDomain );

}

void calculate_pot( struct domain * theDomain ){

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

   double Ptemp = pot;
   for( i=imax ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double dr = c->dr;
      double g = get_g( c ); 
      if( i==imax ) c->pot = Ptemp; else c->pot = Ptemp - .5*g*dr;
      Ptemp -= g*dr;
   }

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

void calculate_fgrav( struct domain * theDomain ){

   int rank=theDomain->rank;
   int size=theDomain->size;

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int i;
   double Int = 0.0;

   int Ng = theDomain->Ng;
   int imax = Nr-Ng;
   if( rank==size-1 ) imax = Nr;
   int imin = Ng;
   if( rank==0 ) imin=0;
   double Itot = 0.0;
   double Itemp = 0.0;

   for( i=imax ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      double dr  = c->dr;
      double rho = c->prim[RHO];
      double v   = c->prim[VRR];
      Itemp -= rho*v*dr;
      c->fgrav = Itemp;
   }

   for( i=imax-1 ; i>=0 ; --i ){
      struct cell * c = theCells+i;
      double dr  = c->dr;
      double rho = c->prim[RHO];
      double v   = c->prim[VRR];
      c->fgrav = Int;
      Int += rho*v*dr;
      if( i==imin ) Itot = Int;
   }

   double Itot_outer = 0.0;
   double Irecv = 0.0;

   int nrk;
   for( nrk=size-1 ; nrk >= 0 ; --nrk ){
      if( rank==nrk ){
         Itot_outer = Irecv;
         Itot += Itot_outer;
         if( nrk > 0 ){
            MPI_Send( &Itot , 1 , MPI_DOUBLE , nrk-1 , 666 , MPI_COMM_WORLD );
         }
      }
      if( rank==nrk-1 ){
         MPI_Status status;
         MPI_Recv( &Irecv , 1 , MPI_DOUBLE , nrk , 666 , MPI_COMM_WORLD , &status );
      }
   }

   for( i=0 ; i<Nr ; ++i ){
      struct cell * c = theCells+i;
      c->fgrav += Itot_outer;
      double M = c->miph;
      double r = c->riph;
      double g = -grav_G*M/r/r;
      c->fgrav *= .5*g;
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
      rm = rp - c->dr;
      double dV = get_dV( rp , rm );
      grav_src( c , dV*dt );
   }

}


