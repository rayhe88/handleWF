#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include "Wavefunction.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <CL/sycl.hpp>

class Field {
public:
  Field(Wavefunction &wf);

  double orbital(Rvector r, int i);
  double gaussiana(Rvector r, int mu);
  double Density(Rvector r);
  void evalDensity();
  void evalDensity2();
  void evalDensity2D();

  void evalDensity_sycl();
  void evalDensity_sycl2();
  static SYCL_EXTERNAL double DensitySYCL(int, int, int*, int*, double*, double*, double *,double*, double*, double*);
  static SYCL_EXTERNAL double DensitySYCL2(int, int, int*, int*, double*, double*, double *,double*, double*, double*);

  void dumpXYZ();

  void dumpCube(double, double, double, double, int, int, int, vector<double>);

private:
  Wavefunction &wf;
};

#endif
