#ifndef MM_MILKYWAY_H
#define MM_MILKYWAY_H
#include "density.h"

class MM_MilkyWay:public density{
 public:
  double getMaxRadius(double z);
  double getDensity(double R, double phi, double theta);

  double get_disc_rho(double RR, double ZZ, double SIGMAD0, double ZD, double RD);
};

#endif
