
#include "../paul.h"

static int NL = 0;
static double * rr  = NULL;
static double * rho = NULL;
static double * Pp  = NULL;
static double * vr  = NULL;
static double * Om  = NULL;

int countlines(char * filename){
   FILE *pFile = fopen(filename, "r");
   int lines=0;
   char c;
   while ((c = fgetc(pFile)) != EOF){
      if (c == '\n') ++lines;
   }
   fclose(pFile);
   return(lines);
}

int getTable( void ){
   int nL = countlines("Initial/initial.dat");
   rr  = (double *) malloc( nL*sizeof(double) );
   rho = (double *) malloc( nL*sizeof(double) );
   Pp  = (double *) malloc( nL*sizeof(double) );
   vr  = (double *) malloc( nL*sizeof(double) );
   Om  = (double *) malloc( nL*sizeof(double) );
   FILE * pFile = fopen("Initial/initial.dat","r");
   int l;
   for( l=0 ; l<nL ; ++l ){
      fscanf(pFile,"%lf %lf %lf %lf %lf\n",&(rr[l]),&(rho[l]),&(Pp[l]),&(vr[l]),&(Om[l]));
   }
   fclose(pFile);
   return(nL);
}

void setICparams( struct domain * theDomain ){
   NL = getTable();
}

void initial( double * prim , double r ){

   int l=0;
   while( rr[l] < r && l < NL-2 ) ++l;
   if( l==0 ) ++l;

   double rp = rr[l];
   double rm = rr[l-1];
   double drm = fabs(r-rm);
   double drp = fabs(rp-r);

   double rh    = (rho[l-1]*drp + rho[l]*drm)/(drp+drm);
   double V     = (vr[l-1]*drp  + vr[l]*drm )/(drp+drm);
   double X = 1.0;

   if( l==NL-2 ){
      V = 0.0;
      X = 0.0;
      rh = 1./r/r;
   }
   double P = rh*1e-8;

   prim[RHO] = rh;
   prim[PPP] = P;
   prim[VRR] = V;
   prim[XXX] = X;
   prim[AAA] = 0.0;

}

void freeTable( void ){
   free(rr);
   free(rho);
   free(Pp);
   free(vr);
   free(Om);
}
