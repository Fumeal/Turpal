#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <utility>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;

double modt (double a, double b, double c, double Mhh)
{
  if (abs(a) <= Mhh)
    return a;

  else if ( (a>=0.) && (b>=0.) && (c>=0.) )
    return min( a , min(b,c) );

  else if ( (a<=0.) && (b<=0.) && (c<=0.) )
    return - min( abs(a) , min(abs(b),abs(c)) );
  
  else
    return 0.;  
}


double modt2 (double a, double b, double Mhh)
{
  if (abs(a) <= Mhh)
    return a;

  else if ( (a>=0.) && (b>=0.) )
    return min( a , b );

  else if ( (a<=0.) && (b<=0.) )
    return - min( abs(a) , abs(b) );
  
  else
    return 0.;  
}




int main ()
{
  //Ordre 4 equation de transport
  double a = 10. ; //vitesse de transport 

  int N = 100; //nombre d'éléments dans la discrétisation.
  vector<vector <double> > U ; //Vecteur contenant notre solution
  U.resize(N);
  double dx = 0.1;
  double dt = 0.99*(1./2.)*dx/a;
  
  double Tmax = 1.;

  //constantes pour la méthode de Galerkine :
  double a0 = 1.;
  double a1 = 12./dx;
  double a2 = 180./(dx*dx);
  double a3 = 14./(23.*dx*dx*dx);
  double I0, I1, I2, I3;
  double hm, hp; // hj+ou-1/2
  double U0i_1;

  double up,um , dpu0i,dmu0i , dmu0ip1,dpu0ip1 , ut,utt , beta;

  //constantes de vérification :
  double TV, norme, CFL;
 
  //conditions limites à gauche :
  double load = 1.;



  //Initialisation du vecteur U:
  for (int i=0; i<(N-N/2);i++){
    U[i].resize(4);
    U[i][0] = 1.;
    U[i][1] = 0.;
    U[i][2] = 0.;
    U[i][3] = 0.;
    cout << i<<endl;
  }  
  for (int i=N/2; i<N;i++){
    U[i].resize(4);
    U[i][0] = 1./2.;
    U[i][1] = 0.;
    U[i][2] = 0.;
    U[i][3] = 0.;
  }

  double Mhh = (2./3.) * (1./2.)/(dx*dx); //à def en fonction de la donnée initiale.

  
  //Faire des tests de type CFL et tout :
  
  double t = 0.;
  while (t<Tmax) //Boucle du schéma
    {
      //def u0_t,u0_tt,...
      //def les trucs pour tester la convergence
      TV = 0.;
      beta = abs(a) ;

      U0i_1 = U[0][0];

      //Calcul du flux à l'interface du bas hm :-----------------
      ut  =   6.*U[0][1] + 30.*U[0][2] + (13./24.)*U[0][3]  ;
      dpu0i = U[1][0] - U[0][0];
      up = U[0][0] - modt2(ut,dpu0i,Mhh);
      
      hm = (1./2.)*( a*load + a*up - beta*(up - load));
      //---------------------------------------------------------
      
      //Calcul du flux à l'interface de droite hp :--------------
      ut   =   6.*U[0][1]   + 30.*U[0][2]   + (13./24.)*U[0][3];
      dpu0i = U[1][0] - U[0][0];
      dmu0i = 2*(U[0][0] - load);
      um = U[0][0] + modt(ut,dpu0i,dmu0i,Mhh);

      utt  = - 6.*U[1][1] + 30.*U[1][2] - (13./24.)*U[1][3];
      dpu0ip1 = U[2][0] - U[1][0];
      dmu0ip1 = U[1][0] - U[0][0];
      up = U[1][0] - modt(utt,dpu0i,dmu0i,Mhh);
      
      hp = (1./2.)*( a*um + a*up - beta*(up - um));
      //---------------------------------------------------------
      
      
      I1 = a*a0*U[0][0]*dx;
      I2 = a*a1*U[0][1]*dx*dx*dx/6.;
      I3 = a*(  a0*U[0][0]*2*dx*dx*dx  +  a2*U[0][2]*pow(dx,5)/5 );
	  
      U[0][0] -= (dt/dx)*(hp - hm);
      U[0][1] += (dt/(dx*dx))*I1 - (dt/dx)*(hp+hm)/2.;
      U[0][2] += (dt/(dx*dx*dx))*I2 - (dt/dx)*(hp-hm)/6.;
      U[0][3] += (dt/(dx*dx*dx*dx))*I3 - (dt/dx)*(hm-hp);
      
      
      for (int i=1;i<N-2;i++)
	{
	  hm=hp;//On réutilise le flux de la maille précédente.

	  //Calcul du flux à l'interface de droite hp :------------------------
	  ut   =   6.*U[i][1]   + 30.*U[i][2]   + (13./24.)*U[i][3];
	  dpu0i = U[i+1][0] - U[i][0];
	  dmu0i = U[i][0] - U0i_1;
	  um = U[i][0] + modt(ut,dpu0i,dmu0i,Mhh);
	  
	  utt  = - 6.*U[i][1] + 30.*U[i][2] - (13./24.)*U[i][3];
	  dpu0ip1 = U[i+2][0] - U[i+1][0];
	  dmu0ip1 = U[i+1][0] - U[i][0];
	  up = U[i+1][0] - modt(utt,dpu0ip1,dmu0ip1,Mhh);
	  
	  hp = (1./2.)*( a*um + a*up - beta*(up - um));
	  //-------------------------------------------------------------------
	  
	  // Calcul de l'integral S ( f(u^h(x,t)) . d/dx(v_l^(j) (x) )dx = I_l pour le cas 
	  I1 = a*a0*U[i][0]*dx;
	  I2 = a*a1*U[i][1]*dx*dx*dx/6.;
	  I3 = a*(  a0*U[i][0]*2*dx*dx*dx  +  a2*U[i][2]*pow(dx,5)/5 );   

	  // Evolution en temps avec un schéma de type Eul_exp
	  U0i_1    =  U[i][0]; //permet de le réutiliser sur la maille suivante.
	  U[i][0] -= (dt/dx)*(hp - hm);
	  U[i][1] += (dt/(dx*dx))*I1 - (dt/dx)*(hp+hm)/2.;
	  U[i][2] += (dt/(dx*dx*dx))*I2 - (dt/dx)*(hp-hm)/6.;
	  U[i][3] += (dt/(dx*dx*dx*dx))*I3 - (dt/dx)*(hm-hp);
	    


	  TV += abs(U[i][0] - U[i-1][0]);

	}
     

      //Calcul des deux derniers termes (condition de bord droit)
      
      //Terme N-2 :-----------------------------------------------------------
      hm = hp;

      ut   =   6.*U[N-2][1] + 30.*U[N-2][2] + (13./24.)*U[N-2][3];
      dpu0i = U[N-1][0] - U[N-2][0];
      dmu0i = U[N-2][0] - U0i_1;
      um = U[N-2][0] + modt(ut,dpu0i,dmu0i,Mhh);

      utt  = - 6.*U[N-1][1] + 30.*U[N-1][2] - (13./24.)*U[N-1][3];
      dmu0ip1 = U[N-1][0] - U[N-2][0];
      up = U[N-1][0] - modt2(utt,dmu0ip1,Mhh);
      
      hp = (1./2.)*( a*um + a*up - beta*(up - um));
      
      I1 = a*a0*U[N-2][0]*dx;
      I2 = a*a1*U[N-2][1]*dx*dx*dx/6.;
      I3 = a*(  a0*U[N-2][0]*2*dx*dx*dx  +  a2*U[N-2][2]*pow(dx,5)/5 );
	  
      U0i_1 = U[N-2][0];
      U[N-2][0] -= (dt/dx)*(hp - hm);
      U[N-2][1] += (dt/(dx*dx))*I1 - (dt/dx)*(hp+hm)/2.;
      U[N-2][2] += (dt/(dx*dx*dx))*I2 - (dt/dx)*(hp-hm)/6.;
      U[N-2][3] += (dt/(dx*dx*dx*dx))*I3 - (dt/dx)*(hm-hp);
      //----------------------------------------------------------------------
      
      
      //Terme N-1 :-----------------------------------------------------------
      hm = hp;
      
      ut   =   6.*U[N-1][1] + 30.*U[N-1][2] + (13./24.)*U[N-1][3];
      dmu0i = U[N-1][0] - U0i_1;
      um = U[N-1][0] + modt2(ut,dmu0i,Mhh);

      hp = a*um; //On a um = up donc hp = a*um

      I1 = a*a0*U[N-1][0]*dx;
      I2 = a*a1*U[N-1][1]*dx*dx*dx/6.;
      I3 = a*(  a0*U[N-1][0]*2*dx*dx*dx  +  a2*U[N-1][2]*pow(dx,5)/5 );
	  
      U0i_1 = U[0][0];
      U[N-1][0] -= (dt/dx)*(hp - hm);
      U[N-1][1] += (dt/(dx*dx))*I1 - (dt/dx)*(hp+hm)/2.;
      U[N-1][2] += (dt/(dx*dx*dx))*I2 - (dt/dx)*(hp-hm)/6.;
      U[N-1][3] += (dt/(dx*dx*dx*dx))*I3 - (dt/dx)*(hm-hp);
      //----------------------------------------------------------------------

      t+=dt;


      //Sauvegarde de la solution : -------------------------------------------
      ofstream flux;
      flux.open("sol"+to_string((int) floor(t/dt))+".dat");

      for (int i=0;i<N;i++)
	{  flux << (i+(1./2.))*dx <<" "<<U[i][0]<<endl ;  }

      flux.close();
      //-----------------------------------------------------------------------
      cout <<t<<endl; 
    }
   


  return 0;
}

