#include "MM_MilkyWay.h"
#include <iostream>
#include <cmath>
using namespace std;

// BULGE
#define solarmass 1.989*pow(10.,30) //kg
#define rhob0 95.6*solarmass*pow(10.,9) //kpc-3
#define q 0.5
#define rcut 2.1 //kpc
#define r0 0.075 //kpc
#define alpha 1.8

// DISCS
#define zd_thick 0.9 //kpc
#define zd_thin 0.3 //kpc
#define rd_thick 3.31 //kpc
#define rd_thin 2.9 //kpc
#define sigmad0_thick 209.5*pow(10.,6)*solarmass
#define sigmad0_thin 816.6*pow(10.,6)*solarmass

//calculated at R=RR and z=ZZ:
double MM_MilkyWay::get_disc_rho(double RR, double ZZ, double SIGMAD0, double ZD, double RD){

  double exp_argument_disc = -(abs(ZZ)/ZD)-(RR/RD);       //protected exponential 
  double disc_rho;
  if (exp_argument_disc>-400){
    disc_rho = (SIGMAD0/2*ZD)*exp(exp_argument_disc);
  }
  else{
    disc_rho = 0.;
  }
  
  return disc_rho;
}


// DARK MATTER HALO
#define rhoh0 0.00846*solarmass*pow(10.,9) //kpc-3  
#define rh 20.2 //kpc







double MM_MilkyWay::getDensity(double R, double phi, double theta){

  double z=R*cos(phi);
  double xx = R/rh;
  double zq2 = (z*z/(q*q));
  double R2 = R*R;

  double distor = sqrt(R2 + zq2);
  double expor = sqrt(R2 + zq2)/rcut;
  double exp_argument_bulge = -(expor*expor);         //protected exponential
  double bulge_rho;
  if(exp_argument_bulge>-400){ 
    bulge_rho=(rhob0/pow(1.0+(distor/r0),alpha))*exp(exp_argument_bulge);
  }
  else{
    bulge_rho=0.;
  }


  

  double thick_disc_rho = get_disc_rho(R,z,sigmad0_thick,zd_thick,rd_thick);
  double thin_disc_rho = get_disc_rho(R,z,sigmad0_thin,zd_thin,rd_thin);

  double discs_rho = thick_disc_rho + thin_disc_rho;

  
  double dmhalo_rho = rhoh0/((xx+0.001)*(1.0+xx)*(1.0 + xx));



////////////////////////////////////////////////////
  
// TOTAL DENSITY DISTRIBUTION
  double rho = bulge_rho + discs_rho + dmhalo_rho;

  return rho;
}    



double MM_MilkyWay::getMaxRadius(double z){
  double maxradius=log((sigmad0_thick)/(2*0.000000001*solarmass*zd_thick*exp(z/zd_thick)))*rd_thick;  //=>density is lower than 0.000000001 solarmasses per kpc-2
  return maxradius;
}
