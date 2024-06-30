#ifndef _PROMOLECULE_HPP_
#define _PROMOLECULE_HPP_

#include "Atom.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <sycl/sycl.hpp>
#include <vector>

class Promolecule {
public:
  Promolecule(std::string fname);
  void loadXYZ(std::string fname);
  void test();

  void dumpCube(vector<double>, std::string filename);

  void dumpData(std::string filename);

  void evalNCI(std::string name);

private:
  int natom;
  int nx, ny, nz;
  double x0, y0, z0;
  double h;
  std::vector<Atom> atoms;
  std::vector<double> den;
  std::vector<double> rdg;
  double density(double x, double y, double z, std::vector<Atom> atoms);
};

#endif
