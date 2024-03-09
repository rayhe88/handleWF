#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include "Wavefunction.hpp"
#include <cmath>
#include <iostream>
#include <vector>

class Field {
public:
  Field(Wavefunction &wf);

  double orbital(Rvector r, int i);
  double gaussiana(Rvector r, int mu);
  double Density(Rvector r);
  void evalDensity();
  void dumpXYZ();
  void dumpCube(double, double, double, double, int, int, int, vector<double>);

private:
  Wavefunction &wf;
};

#endif