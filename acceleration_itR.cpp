#include "acceleration_itR.h"
#include <iostream>
#include <cmath>

#include "MM_MilkyWay.h"

using namespace std;

MM_MilkyWay md;

	void acceleration_itR::getAcc_itR(double r2, double phi2, double theta2,        double r1max, double direction, int N){
	int it=0;
	 
	//double direction  = -1. if radius changes towards centre of galaxy
	//                  = +1. if radius changes towards outside of galaxy

	 //int N = 50;                             //initial point frequency definition. Redefined under phi iteration
	 double dN = ((double)(N));



	 double sphi2 = sin(phi2); double cphi2 = cos(phi2);
	 double stheta2 = sin(theta2);  double ctheta2 = cos(theta2);

	  double R0=r2;  
	  double deltaR = 1.0e-2;//(r1max+R0)/(dN);
	  int itR=1;

	  
	  for (int i=0; i<3; i++){
	    acc[i]=0.;           //initializes acc =0
	  }
	  
	  //varies R
	  double R=R0;
	  double acc_R[3];

	

	double Rcondition;
	if (direction < 0.){
	Rcondition = r2+1e-15;
	}
	if (direction > 0.){
	Rcondition = abs(r2+r1max);
	}

	  while ((R < Rcondition) && (R>0)) {
	      for (int i=0; i<3; i++){         
		acc_R[i]=0.;           //initializes fixed R acc =0
	      }
	      
	    //varies phi
	      for (int j_phi=0; j_phi<=N; j_phi++){
		double phi=((double)(j_phi))*pi/dN;
		double acc_R_phi[3];
		for (int i=0; i<3; i++){
		  acc_R_phi[i]=0.;       //initializes fixed R and phi acc =0
		}
		double sphi = sin(phi); double cphi = cos(phi);


		double z1=R*cphi+r2*cphi2;                  //this thing makes the theta density higher when close to "equator"
		  int NN=0;
		  double dNN=0.;
		  double ddNN=0.;

		/*if(abs(z1)<15.){
		ddNN = (2./45.)*z1*z1 + (-4./3.)*z1 + 10.*N; 
		NN = ((int)(ddNN)); 
		dNN = ((double)(NN));
		}
		if(abs(z1)>15.){
		NN = 0;
		dNN = 0.0;
		}*/



		  if(abs(z1)<=1.0){
		    NN=10*N;
		    dNN=((double)(NN));
		  }
	      	  if((abs(z1)>1.0) && (abs(z1)<=2.0)){
		    NN=2*N;
		    dNN=((double)(NN));
		  }
		  
		  if(abs(z1)>2.0){
		    NN=N;
		    dNN=((double)(NN));
		  }
		  if(abs(z1)>=10.0){
		    NN=N/10;
		    dNN=((double)(NN));
		  }


                  if(j_phi==0 || j_phi==N){
		    NN=1;
		    dNN=1.;
		  }


	
		//varies theta
		for (int j_theta=0; j_theta<=NN-1; j_theta++){
		  double theta=theta2+((2.0*((double)(j_theta))+1.0)/2.0)*pi/dNN;       //avoids divergence at origin
		  double stheta = sin(theta); double ctheta = cos(theta);

		  //ref frame change to find density
		  double x1=R*sphi*ctheta+r2*sphi2*ctheta2;
		  double y1=R*sphi*stheta+r2*sphi2*stheta2;
		  //double z1=R*cphi+r2*cphi2;

		  double x1_opposed=R*sphi*cos(theta+pi)+r2*sphi2*ctheta2;
		  double y1_opposed=R*sphi*sin(theta+pi)+r2*sphi2*stheta2;
		  double z1_opposed=R*cphi+r2*cphi2;

		  double r1=sqrt(x1*x1+y1*y1+z1*z1);
	      	  double phi1=acos(z1/r1);
		  double theta1=atan2(y1,x1);

		  double r1_opposed=sqrt(x1_opposed*x1_opposed+y1_opposed*y1_opposed+z1_opposed*z1_opposed);
		  double phi1_opposed=acos(z1_opposed/r1_opposed);
		  double theta1_opposed=atan2(y1_opposed,x1_opposed);
		  

		  it++;

		    double pointDensity=md.getDensity(r1, phi1, theta1)-md.getDensity(r1_opposed, phi1_opposed, theta1_opposed);

		  
		  double directionVector[3];
		  directionVector[0] = sphi*ctheta;       //for acc to be found in cartesian coordinates
		  directionVector[1] = sphi*stheta;       
		  directionVector[2] = cphi;

		  double acc_point[3];
		  for (int i=0; i<3; i++){
		    acc_point[i]=0.;       //initializes fixed R and phi and theta acc =0
		  }	  
		  for (int i=0; i<3; i++){
		    acc_point[i]=G*pointDensity*directionVector[i];
		    acc_R_phi[i]+=acc_point[i]*(pi/dNN)*sphi;    //adds point term to fixed R and phi acc
		  }
		}  // end for theta

		for (int i=0; i<3; i++){
		  acc_R[i]+=acc_R_phi[i]*pi/(dN+1.0);          //adds circumference term to fixed R acc
		}
	      }  // end for phi


	///////////////// ADAPTIVE RADIUS //////////////////////////////////
	      double accd = 0.0; double accdel = 0.0;

	      for (int i = 0; i < 3; i++) {
		accd += acc[i]*acc[i]; accdel += acc_R[i]*acc_R[i];
	      }
	      accdel *= deltaR;
	      if (accdel < accd*1.0e-6) {
		deltaR*=1.4;                                //increases deltaR for improved efficiency
	//      	cout<<R<<"  increasing\n";
	      }                    
	      if ((itR > 2) && (accdel > accd*1.0e-4)) {
		if (deltaR > 0.001) {
	      	deltaR*=0.7;                                //decreases deltaR for improved accuracy
	//	cout<<R<<"  decreasing " << deltaR << " "<< accdel << " " << accd << "\n";
		}
	////////////////////////////////////////////////////////////////////									     
	      }


		for (int i=0; i<3; i++){
		  acc[i]+=acc_R[i]*deltaR;                  //adds sphere term to total acc
		}
		itR++;


		//	cout<<"  R= "<<R << " " << accdel/accd <<endl;
	      R += direction * deltaR;

	  }

//cout<<"\n   acc in "<<direction<<" direction is     "<<acc[0]<<"  "<<acc[1]<<"  "<<acc[2]<<endl<<endl;
	}
