
#include "paul.h"

void report( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   struct cell * theCells = theDomain->theCells;
   int rank = theDomain->rank;
   int size = theDomain->size;
   double t = theDomain->t;

   double Mach_Avg = 0.0;
   double Mass   = 0.0;
   double r_dens = 0.0;
   int imin = Ng;
   int imax = Nr-Ng;
   if( rank==0 ) imin = 0;
   if( rank==size-1 ) imax = Nr;

   int i;
   for( i=imin ; i<imax ; ++i ){
      double rho = theCells[i].prim[RHO];
      double r   = theCells[i].riph;
      if( rho > 0.5 ) r_dens = r;
      double v   = theCells[i].prim[VRR];
      double P   = theCells[i].prim[PPP];
      double cs  = sqrt(1.4*P/rho);
      Mass += theCells[i].cons[DDD];
      Mach_Avg += fabs(v)/cs*theCells[i].cons[DDD];
      
   }

   MPI_Allreduce( MPI_IN_PLACE , &r_dens , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &Mass , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &Mach_Avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );

   Mach_Avg /= Mass;

   if( rank==0 ){
      FILE * pFile = fopen("report.dat","a");
      fprintf(pFile,"%e %e %e\n", t , r_dens , Mach_Avg );
      fclose(pFile);
   }

}


