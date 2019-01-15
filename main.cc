#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <chrono>

using namespace std;

int main ()
{
  //Ordre 4 equation de transport
  double a = 12. ; //vitesse de transport 

  int N = 100; //nombre d'éléments dans la discrétisation.
  vector<vector <double> > U ; //Vecteur contenant notre solution
  U.resize(N);
  double dt = 0.01;
  double dx = 0.01;
  double Tmax = 10;

  //constantes pour la méthode de Galerkine :
  double a0 = 1.;
  double a1 = 12./dx;
  double a2 = 180./(dx*dx);
  double a3 = 14./(23.*dx*dx*dx);
  double I0, I1, I2, I3;
  double hm, hp; // hj+ou-1/2  


  //Initialisation du vecteur U:
  for (int i=0; i<N;i++){
    U[i].resize(4);
    U[0] = 0.;
    U[1] = 0.;
    U[2] = 0.;
    U[3] = 0.;
  }

  //Faire des tests de type CFL et tout :
  
  double t = 0.;
  while (t<Tmax) //Boucle du schéma
    {
      //def u0_t,u0_tt,...
      //def les trucs pour tester la convergence

      hm = 0;//à def pour les conditions au bord. 
     
      for (int i=0;i<N-1;i++)
	{
	  hm = hp;
	  //Faire le calcul de hp


	  // Calcul de l'integral S ( f(u^h(x,t)) . d/dx(v_l^(j) (x) )dx = I_l pour le cas 
	  I1 = a*a0*U[i][0]*dx;
	  I2 = a*a1*U[i][1]*dx*dx*dx/6.;
	  I3 = a*(  a0*U[i][0]*2*dx*dx*dx  +  a2*U[i][2]*pow(dx,5)/5 );   


	  U[i][0] -= (dt/dx)*(hp - hm);
	  U[i][1] += (dt/(dx*dx))*I1 - (dt/dx)*(hp+hm)/2.;
	  U[i][2] += (dt/(dx*dx*dx))*I2 - (dt/dx)*(hp-hm)/6.;
	  U[i][3] += (dt/(dx*dx*dx*dx))*I3 - (dt/dx)*(hm-hp);
	    


	}
     

      //Calcul dernier terme (condition de bord et tout



      t+=dt;
    }
   






  return 0;
}
