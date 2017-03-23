
#include "paul.h"

void calculate_pot( struct domain * );

void report( struct domain * theDomain ){

   calculate_pot( theDomain );

   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   struct cell * theCells = theDomain->theCells;
   int rank = theDomain->rank;
   int size = theDomain->size;
   double t = theDomain->t;

   double Mach_Avg = 0.0;
   double Mass   = 0.0;
   double XMass   = 0.0;
   double Turb = 0.0;
   double r_dens = 0.0;
   double r_shock = 0.0;
   double E_Hydro = 0.0;
   double E_Kin   = 0.0;
   double E_Grav  = 0.0;
   int imin = Ng;
   int imax = Nr-Ng;
   if( rank==0 ) imin = 0;
   if( rank==size-1 ) imax = Nr;

   double w_avg   = 0.0;
   double rho_avg = 0.0;
   double P_avg = 0.0;
   double r_avg   = 0.0;
   double r2_avg  = 0.0;
   double v_avg   = 0.0;
   double a_avg   = 0.0;
   double grho_avg= 0.0;
   double gP_avg  = 0.0;

   int i;
   for( i=imin ; i<imax ; ++i ){
      double rho = theCells[i].prim[RHO];
      double r   = theCells[i].riph;
      double dr  = theCells[i].dr;

      if( rho > 0.5 ) r_dens = r;
      double v   = theCells[i].prim[VRR];
      double P   = theCells[i].prim[PPP];
      double cs  = sqrt(1.4*P/rho);
      if( fabs(v) > cs ) r_shock = r;
      Mass += theCells[i].cons[DDD];
      XMass += theCells[i].cons[XXX];
      Turb += theCells[i].cons[AAA];
      Mach_Avg += fabs(v)/cs*theCells[i].cons[DDD];
      E_Hydro += theCells[i].cons[TAU];
      E_Kin += .5*theCells[i].cons[DDD]*v*v;
      E_Grav  += .5*theCells[i].cons[DDD]*theCells[i].pot;

      double alpha = theCells[i].prim[AAA];
      double X   = theCells[i].prim[XXX];
      double w = X*(1.-X);
      
      double grho = theCells[i].grad[RHO];
      double gP   = theCells[i].grad[PPP];

      w_avg   += w*dr;
      rho_avg += rho*w*dr;
      P_avg   += P*w*dr;
      r_avg   += r*w*dr;
      r2_avg  += r*r*w*dr;
      v_avg   += v*w*dr;
      a_avg   += alpha*w*dr;
      grho_avg += grho*w*dr;
      gP_avg  += gP*w*dr;

   }

   MPI_Allreduce( MPI_IN_PLACE , &r_dens , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &r_shock , 1 , MPI_DOUBLE , MPI_MAX , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &Mass , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &XMass , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &Turb , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &E_Hydro , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &E_Kin , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &E_Grav , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &Mach_Avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );


   MPI_Allreduce( MPI_IN_PLACE , &w_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &rho_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &P_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &r_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &r2_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &v_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &a_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &grho_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
   MPI_Allreduce( MPI_IN_PLACE , &gP_avg , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );

   rho_avg /= w_avg;
   P_avg /= w_avg;
   r_avg /= w_avg;
   r2_avg /= w_avg;
   v_avg /= w_avg;
   a_avg /= w_avg;
   grho_avg /= w_avg;
   gP_avg /= w_avg;

   Mach_Avg /= Mass;

   if( rank==0 ){
      FILE * pFile = fopen("report.dat","a");
      fprintf(pFile,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", t , r_dens , Mach_Avg , E_Hydro , E_Grav , r_avg , r2_avg , rho_avg , v_avg , a_avg , grho_avg , gP_avg , P_avg , Mass , XMass , E_Kin , Turb , r_shock );
      fclose(pFile);
   }

}


