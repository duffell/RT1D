
#include "paul.h"

int mpiSetup( struct domain * theDomain , int argc, char * argv[] ){

   MPI_Comm_size(MPI_COMM_WORLD,&(theDomain->size));
   MPI_Comm_rank(MPI_COMM_WORLD,&(theDomain->rank));

   return(0);

}

