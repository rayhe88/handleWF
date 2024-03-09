#include "Field.hpp"
#include "Atom.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

Field::Field(Wavefunction &wf) : wf(wf) {}

double Field::gaussiana(Rvector r, int mu) {
  int lx = wf.vang[3 * mu];
  int ly = wf.vang[3 * mu + 1];
  int lz = wf.vang[3 * mu + 2];
  Rvector R(wf.atoms[wf.icntrs[mu]].getCoors());

  double x_part = r.get_x() - R.get_x();
  double y_part = r.get_y() - R.get_y();
  double z_part = r.get_z() - R.get_z();

  double diff2 = pow(x_part, 2) + pow(y_part, 2) + pow(z_part, 2);

  double gauss = pow(x_part, lx) * pow(y_part, ly) * pow(z_part, lz);

  gauss *= exp(-wf.depris[mu] * diff2);

  return gauss;
}

double Field::orbital(Rvector r, int i) {
  double orb;

  orb = 0.;
  for (int mu = 0; mu < wf.npri; mu++)
    orb += (wf.dcoefs[i * wf.npri + mu] * gaussiana(r, mu));

  return orb;
}
double Field::Density(Rvector r) {
  double den;
  den = 0;
  for (int i = 0; i < wf.norb; i++)
    den += wf.dnoccs[i] * orbital(r, i) * orbital(r, i);

  return den;
}

void Field::evalDensity() {

  vector<double> field;
  double xmin = -10.0, xmax = 10.0;
  double ymin = -10.0, ymax = 10.0;
  double zmin = -5.0, zmax = 5.0;
  double delta = 0.05;

  int npoints_x = int((xmax - xmin) / delta);
  int npoints_y = int((ymax - ymin) / delta);
  int npoints_z = int((zmax - zmin) / delta);

  for (int i = 0; i < npoints_x; i++) {
    double x = xmin + i * delta;
    for (int j = 0; j < npoints_y; j++) {
      double y = ymin + j * delta;
      for (int k = 0; k < npoints_z; k++) {
        double z = zmin + k * delta;
        Rvector r(x, y, z);
        double den = Density(r);

        field.push_back(den);
      }
    }
  }

  dumpCube(xmin, ymin, zmin, delta, npoints_x, npoints_y, npoints_z, field);
  dumpXYZ();
}

void Field::dumpCube(double xmin, double ymin, double zmin, double delta,
                     int nx, int ny, int nz, vector<double> field) {
  std::ofstream fout("density.cube");
  if (fout.is_open()) {

    fout << "Density" << std::endl;
    fout << "By handleWF project" << std::endl;
    fout << std::setw(5) << std::fixed << wf.natm << std::setw(13)
         << std::setprecision(6) << std::fixed << xmin << ymin << zmin
         << std::endl;

    fout << std::setw(5) << std::fixed << nx;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << delta << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0;
    fout << std::endl;

    fout << std::setw(5) << std::fixed << ny;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << delta << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0;
    fout << std::endl;

    fout << std::setw(5) << std::fixed << nz;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << delta;
    fout << std::endl;
  }

  for (auto atom : wf.atoms) {
    fout << std::setw(5) << std::fixed << atom.get_atnum();
    fout << std::setw(13) << std::setprecision(6) << std::fixed
         << atom.get_charge() << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << atom.get_x()
         << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << atom.get_y()
         << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << atom.get_z();
    fout << std::endl;
  }

  int cnt = 0;
  for (auto valor : field) {
    cnt++;
    fout << std::setw(13) << std::setprecision(6) << std::fixed
         << std::scientific << valor;
    if (cnt == 6) {
      fout << std::endl;
    }
  }
  if (cnt != 0)
    fout << std::endl;

  fout.close();
}

void Field::dumpXYZ() {
  std::ofstream fout("structura.xyz");
  if (fout.is_open()) {

    fout << std::setw(4) << std::fixed << wf.natm << std::endl;
    fout << " File created by handleWF code" << std::endl;
    for (auto atom : wf.atoms) {
      fout << std::setw(4) << std::fixed << atom.getSymbol();
      fout << std::setw(13) << std::setprecision(6) << std::fixed
           << atom.get_x() * 0.529177 << ' ';
      fout << std::setw(13) << std::setprecision(6) << std::fixed
           << atom.get_y() * 0.529177 << ' ';
      fout << std::setw(13) << std::setprecision(6) << std::fixed
           << atom.get_z() * 0.529177;
      fout << std::endl;
    }
    fout.close();
  }
}