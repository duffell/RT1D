
#include "../paul.h"

double get_moment_arm( double , double );
double get_dV( double , double );

void output( struct domain * theDomain , char * filestart ){

   struct cell * theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Ng = theDomain->Ng;
   int rank = theDomain->rank;
   int size = theDomain->size;

   char filename[256];
   sprintf(filename,"%s.dat",filestart);

   if( rank==0 ){
      FILE * pFile = fopen( filename , "w" );
      fprintf(pFile,"#r           dr           Density      Pressure     Velocity     X            Alpha\n");
      fclose(pFile);
   }
   MPI_Barrier( MPI_COMM_WORLD );

   int i_min = 0;
   int i_max = Nr;

   if( rank != 0      ) i_min = Ng;
   if( rank != size-1 ) i_max = Nr-Ng;

   int rk;
   for( rk=0 ; rk<size ; ++rk ){
      if( rank==rk ){
         FILE * pFile = fopen( filename , "a" );
         int i,q;
         for( i=i_min ; i<i_max ; ++i ){
            struct cell * c = theCells+i;
            double rp = c->riph;
            double dr = c->dr; 
            double rm = rp-dr;
            double r  = get_moment_arm( rp , rm );
            fprintf(pFile,"%e %e ",r,dr);
            for( q=0 ; q<NUM_Q ; ++q ){
               fprintf(pFile,"%e ",c->prim[q]);
            }
            fprintf(pFile,"%e ",c->miph);
            fprintf(pFile,"\n");
         }
         fclose( pFile );
      }
      MPI_Barrier( MPI_COMM_WORLD );
   }

}
