#ifndef ACCELERATION_ITR
#define ACCELERATION ITR

#define pi 3.14159265359
#define G 6.6740831e-11

class acceleration_itR{
 public:

  double acc[3];
  void getAcc_itR(double r2, double phi2, double theta2,       double r1max, double direction, int N);

};

#endif
