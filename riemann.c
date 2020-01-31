enum{_HLL_,_HLLC_};

#include "paul.h"

static int riemann_solver = 0;
static int rt_flag = 0;
static double gamma_law = 1.0;
static int grav_e_mode = 0;
static int Bal = 0;
static double grav_G = 1.0;

void setRiemannParams( struct domain * theDomain ){
   riemann_solver = theDomain->theParList.Riemann_Solver;
   rt_flag = theDomain->theParList.rt_flag;
   gamma_law = theDomain->theParList.Adiabatic_Index;
   grav_e_mode = theDomain->theParList.grav_e_mode;
   Bal = theDomain->theParList.grav_bal;
   grav_G = theDomain->theParList.grav_G;
}

void prim2cons( double * , double * , double , double );
void flux( double * , double * );
void getUstar( double * , double * , double , double );
void vel( double * , double * , double * , double * , double * , double );
double get_eta( double * , double * , double );

double get_g( struct cell * );
double get_GMr( struct cell * );

void riemann( struct cell * cL , struct cell * cR, double r , double dAdt ){

   double primL[NUM_Q];
   double primR[NUM_Q];

   double drL = .5*cL->dr;
   double drR = .5*cR->dr;

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + cL->grad[q]*drL;
      primR[q] = cR->prim[q] - cR->grad[q]*drR;
   }

   if( Bal ){
      //Add HSE pressure gradient to pL and pR.
      double gHSE_L = cL->prim[RHO]*get_g(cL);
      double gHSE_R = cR->prim[RHO]*get_g(cR);

      primL[PPP] += gHSE_L*drL;
      primR[PPP] -= gHSE_R*drR;
   }

   double Sl,Sr,Ss;

   vel( primL , primR , &Sl , &Sr , &Ss , r );


   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];

   double Flux[NUM_Q];

   double w = cL->wiph;

   if( w < Sl ){
      flux( primL , Fl );
      prim2cons( primL , Ul , 0.0 , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fl[q] - w*Ul[q];
      }
   }else if( w > Sr ){
      flux( primR , Fr );
      prim2cons( primR , Ur , 0.0 , 1.0 );

      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] = Fr[q] - w*Ur[q];
      }
   }else{
      if( riemann_solver == _HLL_ ){
         double Fstar;
         double Ustar;
         double aL =  Sr;
         double aR = -Sl;
 
         prim2cons( primL , Ul , 0.0 , 1.0 );
         prim2cons( primR , Ur , 0.0 , 1.0 );
         flux( primL , Fl );
         flux( primR , Fr );

         for( q=0 ; q<NUM_Q ; ++q ){
            Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
            Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );

            Flux[q] = Fstar - w*Ustar;
         }
      }else{
         double Ustar[NUM_Q];
         double Uk[NUM_Q];
         double Fk[NUM_Q];
         if( w < Ss ){
            prim2cons( primL , Uk , 0.0 , 1.0 );
            getUstar( primL , Ustar , Sl , Ss ); 
            flux( primL , Fk ); 

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
            }    
         }else{
            prim2cons( primR , Uk , 0.0 , 1.0 );
            getUstar( primR , Ustar , Sr , Ss ); 
            flux( primR , Fk ); 

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
            } 
         } 
      }
   } 
 
   if( grav_e_mode == 3 ){
      double G = grav_G;
      double m = cL->miph;
      double r = cL->riph;
      Flux[TAU] += -G*m/r*Flux[RHO];
   }

   if( rt_flag ){
      double prim[NUM_Q];
      double consL[NUM_Q];
      double consR[NUM_Q];
      prim2cons( cL->prim , consL , 0.0 , 1.0 );
      prim2cons( cR->prim , consR , 0.0 , 1.0 );
      double gprim[NUM_Q];
      double gcons[NUM_Q];
      for( q=0 ; q<NUM_Q ; ++q ){
         prim[q] = .5*(primL[q]+primR[q]);
         gprim[q] = (cR->prim[q] - cL->prim[q])/(drL+drR);
         gcons[q] = (consR[q] - consL[q])/(drL+drR);
      }
//If new model, gcons[q] = P^1/gamma*grad( cons/P^1/gamma ).
/*
      double Pgam  = pow( prim[PPP] , 1./gamma_law );
      double PgamL = pow( cL->prim[PPP] , 1./gamma_law );
      double PgamR = pow( cR->prim[PPP] , 1./gamma_law );
      for( q=0 ; q<NUM_Q ; ++q ){
         gcons[q] = Pgam*( consR[q]/PgamR - consL[q]/PgamL )/(drL+drR); 
      }
*/
////////
      double eta = get_eta( prim , gprim , r );
      for( q=0 ; q<NUM_Q ; ++q ){
         Flux[q] += -eta*gcons[q];
      }
   }

   for( q=0 ; q<NUM_Q ; ++q ){
      cL->cons[q] -= Flux[q]*dAdt;
      cR->cons[q] += Flux[q]*dAdt;
   }

}


