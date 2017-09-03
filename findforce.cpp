#include <iostream>
#include <cmath>
#include <string>

#include "density.h"
#include "PointMass.h"
#include "MM_MilkyWay.h"

#include "acceleration_itR.h"

using namespace std;

#define pi 3.14159265359
#define G 6.6740831e-11

MM_MilkyWay MD;

acceleration_itR inside;
acceleration_itR outside;


void getAcc(double r2, double phi2, double theta2,        double r1max,    int N){

inside.getAcc_itR(r2, phi2, theta2, r1max, -1., N);
//cout<<"\n**************changing directions**************\n";
outside.getAcc_itR(r2, phi2, theta2, r1max, +1., N);

double acc[3];
for (int i = 0; i < 3; i++) {
acc[i] = inside.acc[i] + outside.acc[i];
}

  cout<<"r2="<<r2<<"  phi2="<<phi2<<"  theta2="<<theta2<<"   acceleration(x,y,z)_is:  "<<acc[0]<<"  "<< acc[1]<< "  "<< acc[2]<<endl;

}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



int main(){

/*double x = 15.; double y = 0; 
double rsqxy = x*x + y*y; double theta=atan2(y,x);
for (double z = 0.; z <= 100.; z += 0.1){
	  double rsq = rsqxy + z*z;
	  double R = sqrt(rsq);
      	  double phi = acos(z/R);
double phi=0.5*pi;
double sphi=sin(phi); double cphi=cos(phi);
double theta=0;
double stheta=sin(theta); double ctheta=cos(theta);    //(this is all messed up)
for(double R=1.; R<=500; R+=1){
     cout<<z<<"    "<< MD.getDensity(R,phi,theta)<<endl;
}*/



for(int N = 10; N < 500; N += 10){
cout<<"N= "<<N<<"   ";
getAcc(8, pi*0.5, 0, 30, N);
cout<<endl<<endl;
}




//for(double r1max=25.; r1max<=200.; r1max+=5.){
//  //  double r1max = 25.;
//    cout<<"r1max= "<<r1max<<"     ";
//  getAcc(8, pi/2., 0,    r1max, 50);
//  cout<<"\n\n";
// }



/*double r1max = 20.; 
double r2 = 8;
double phi2 = pi/2.;
double sphi2=sin(phi2);
for(double n=0.; n<=2.; n+=0.1){ 
 cout<<"x2= "<<r2*sphi2*cos(n*pi)<<"  y2= "<< r2*sphi2*sin(n*pi)<<" ||| ";
 getAcc(r2, phi2, n*pi, r1max);
}*/

  

}
