#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <cmath>
#include <iostream>
#include <string>
#include <sycl/sycl.hpp>
#include <vector>

#include "Wavefunction.hpp"

class Field {
   public:
    Field(Wavefunction &wf);

    double orbital(std::vector<double> r, int i);
    double gaussiana(std::vector<double> r, int mu);
    double Density(std::vector<double> r);
    void evalDensity();
    void evalDensity2();
    void evalDensity2D();

    void evalDensity_sycl();
    void evalDensity_sycl2();
    static SYCL_EXTERNAL double DensitySYCL(int, int, int *, int *, double *,
                                            double *, double *, double *,
                                            double *, double *);
    static SYCL_EXTERNAL double DensitySYCL2(int, int, const int *, const int *,
                                             const double *, const double *,
                                             const double *, const double *,
                                             const double *);

    void dumpXYZ(std::string filename);

    void dumpCube(double, double, double, double, int, int, int, vector<double>,
                  std::string filename);

   private:
    Wavefunction &wf;
};

#endif
