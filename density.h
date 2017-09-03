#ifndef DENSITY_H
#define DENSITY_H

class density{
 public:
  double rho;
  virtual double getDensity(double R, double phi, double theta)=0;

};
#endif
