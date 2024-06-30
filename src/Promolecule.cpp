#include "Promolecule.hpp"
#include "Atom.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sycl/sycl.hpp>

inline bool is_number(const string &s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

Promolecule::Promolecule(std::string fname) {
  h = 0.1;
  loadXYZ(fname);
}

void Promolecule::loadXYZ(std::string fname) {
  int natm;
  char tmp[10];
  string s;
  double xt, yt, zt;

  std::vector<double> coorx;
  std::vector<double> coory;
  std::vector<double> coorz;

  string line;
  ifstream finp;

  finp.open(fname.c_str(), std::ios::in);
  if (!finp.good()) {
    cout << " The file [" << fname << "] can't be opened!" << endl;
  }
  finp.seekg(finp.beg);

  getline(finp, line); // read the line with the number of nuclei
  sscanf(line.c_str(), "%d", &natm);
  getline(finp, line); // read the comment line;
  for (int i = 0; i < natm; i++) {
    getline(finp, line); // read the line with information of centres
    sscanf(line.c_str(), "%s %lf %lf %lf", tmp, &xt, &yt, &zt);

    xt *= 1.8897259886;
    yt *= 1.8897259886;
    zt *= 1.8897259886;

    coorx.push_back(xt);
    coory.push_back(yt);
    coorz.push_back(zt);

    std::sort(coorx.begin(), coorx.end());
    std::sort(coory.begin(), coory.end());
    std::sort(coorz.begin(), coorz.end());
    double xmin = coorx.front() - 3.;
    double xmax = coorx.back() + 3.;
    double ymin = coory.front() - 3.;
    double ymax = coory.back() + 3.;
    double zmin = coorz.front() - 3.;
    double zmax = coorz.back() + 3.;

    nx = static_cast<int>(std::ceil((xmax - xmin) / h));
    ny = static_cast<int>(std::ceil((ymax - ymin) / h));
    nz = static_cast<int>(std::ceil((zmax - zmin) / h));
    x0 = xmin;
    y0 = ymin;
    z0 = zmin;

    if (!is_number(tmp)) {
      s = string(tmp);
      s.erase(std::remove_if(s.begin(), s.end(),
                             [](char ch) { return std::isdigit(ch); }),
              s.end());
      atoms.push_back(Atom(s, xt, yt, zt));
    } else {
      atoms.push_back(Atom(atoi(tmp), xt, yt, zt));
    }
  }
  natom = atoms.size();

  finp.close();
}

void Promolecule::test() {
  std::cout << " Test of load XYZ" << std::endl;
  for (auto atom : atoms) {
    std::cout << setw(10) << atom.getSymbol();
    std::cout << setw(10) << fixed << setprecision(6) << atom.get_x();
    std::cout << setw(10) << fixed << setprecision(6) << atom.get_y();
    std::cout << setw(10) << fixed << setprecision(6) << atom.get_z()
              << std::endl;
  }
}

double Promolecule::density(double x, double y, double z,
                            std::vector<Atom> atoms) {
  std::vector<double> c1 = {
      0.0000000E+00, 4.6568633E-01, 1.9613583E+01, 1.0682467E+02, 4.9307260E-01,
      2.5663419E+02, 4.0348178E+02, 6.3406713E+02, 9.5831026E+02, 1.6273443E+01,
      1.9788394E+03, 3.5289679E+03, 1.0046303E+02, 6.2991738E+03, 8.4532518E+03,
      1.1107832E+04, 1.4348550E+04, 6.0301596E+02, 8.0816108E+02};
  std::vector<double> c2 = {
      0.0000000E+00, 4.6568622E-01, 1.9568763E+01, 1.0681231E+02, 5.2439107E+00,
      1.1256576E+00, 2.2620806E+00, 6.3404543E+02, 8.9583281E+00, 1.3870221E+03,
      1.9614363E+03, 3.5302178E+03, 4.5960222E+03, 1.5274864E+02, 8.4342635E+03,
      3.1824165E+02, 4.4165626E+02, 1.8325620E+04, 2.3005524E+04};
  std::vector<double> c3 = {
      0.0000000E+00, 3.0019989E+00, 1.3983775E+00, 3.3342545E+01, 3.1734662E+02,
      2.5867085E+02, 4.0837497E+02, 4.6765606E+00, 9.5832030E+02, 1.4029674E+03,
      2.8279579E+01, 7.4670217E+01, 4.6015044E+03, 6.2957609E+03, 2.2368762E+02,
      1.1103799E+04, 1.4364194E+04, 1.8305210E+04, 2.3084071E+04};
  std::vector<double> a1 = {
      0.0000000E+00, 1.9888088E+00, 3.8987496E+00, 6.2028021E+00, 7.2657186E-01,
      9.7173881E+00, 1.1698438E+01, 1.3677050E+01, 1.5655698E+01, 1.5948420E+00,
      1.9625025E+01, 2.2086110E+01, 2.6605426E+00, 2.6057308E+01, 2.8120183E+01,
      3.0186892E+01, 3.2254026E+01, 4.2979610E+00, 4.6484753E+00};
  std::vector<double> a2 = {
      0.0000000E+00, 1.9888088E+00, 3.8982540E+00, 6.2023193E+00, 4.1739221E+00,
      8.3058204E-01, 9.6295177E-01, 1.3677045E+01, 1.3642293E+00, 1.7635609E+01,
      1.9625014E+01, 2.2085863E+01, 2.3998338E+01, 2.9772425E+00, 2.8120238E+01,
      3.6228839E+00, 3.9540388E+00, 3.4332866E+01, 3.6415201E+01};
  std::vector<double> a3 = {
      0.0000000E+00, 1.9888228E+00, 1.6012664E+00, 3.8061318E+00, 7.8497277E+00,
      9.7173599E+00, 1.1698438E+01, 1.1538421E+00, 1.5655807E+01, 1.7635595E+01,
      1.8462843E+00, 2.4886290E+00, 2.3998352E+01, 2.6057363E+01, 3.2971353E+00,
      3.0186664E+01, 3.2254046E+01, 3.4332897E+01, 3.6415230E+01};

  int iatm;
  double xatm, yatm, zatm;
  double r;

  double rho = 0.;
  for (auto atom : atoms) {
    iatm = atom.get_atnum();
    xatm = atom.get_x();
    yatm = atom.get_y();
    zatm = atom.get_z();
    r = sqrt((x - xatm) * (x - xatm) + (y - yatm) * (y - yatm) +
             (z - zatm) * (z - zatm));

    rho += c1[iatm] * exp(-a1[iatm] * r) + c2[iatm] * exp(-a2[iatm] * r) +
           c3[iatm] * exp(-a3[iatm] * r);
  }

  return rho;
}

void Promolecule::evalNCI(std::string name) {
  double x, y, z;
  double rho;
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        x = x0 + i * h;
        y = y0 + j * h;
        z = z0 + k * h;
        rho = density(x, y, z, atoms);
        den.push_back(rho);
      }
    }
  }
  dumpCube(den, name);
}

void Promolecule::dumpCube(std::vector<double> field, std::string filename) {
  std::ofstream fout(filename);
  if (fout.is_open()) {

    fout << "Density" << std::endl;
    fout << "By handleWF project" << std::endl;
    fout << std::setw(5) << std::fixed << natom;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << x0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << y0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << z0;
    fout << std::endl;

    fout << std::setw(5) << std::fixed << nx;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << h << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0;
    fout << std::endl;

    fout << std::setw(5) << std::fixed << ny;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << h << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0;
    fout << std::endl;

    fout << std::setw(5) << std::fixed << nz;
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << 0.0 << ' ';
    fout << std::setw(13) << std::setprecision(6) << std::fixed << h;
    fout << std::endl;
  }

  for (auto atom : atoms) {
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
    fout << std::setw(15) << std::setprecision(6) << std::fixed
         << std::scientific << valor;
    if (cnt == 6) {
      fout << std::endl;
      cnt = 0;
    }
  }
  if (cnt != 0)
    fout << std::endl;

  fout.close();
}