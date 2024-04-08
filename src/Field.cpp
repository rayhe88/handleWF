#include "Field.hpp"
#include "Atom.hpp"
#include "Wavefunction.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <CL/sycl.hpp>

Field::Field(Wavefunction &wf) : wf(wf) {}

double Field::DensitySYCL(int norb, int npri,
                          int *icnt, int *vang,
                          double *r, double *coor, double *depris,
                          double *nocc, double *coef, double *moi){
  double den;
  double x = r[0];
  double y = r[1];
  double z = r[2];

  for(int j=0; j < npri; j++){
    const int vj = 3*j;
    const int centroj = 3*icnt[j];
    const double difx = x - coor[centroj];
    const double dify = y - coor[centroj+1];
    const double difz = z - coor[centroj+2];
    const double rr = difx*difx + dify*dify + difz*difz;

    const double expo = exp(-depris[j]*rr);
    const double lx   = vang[vj];
    const double ly   = vang[vj+1];
    const double lz   = vang[vj+2];
    const double facx = pow(difx,lx);
    const double facy = pow(dify,ly);
    const double facz = pow(difz,lz);

    moi[j]   = facx*facy*facz*expo;
  }


  den=0.0;
  for(int i=0; i < norb; i++){
     double mo=0.0;
     const int i_prim= i*npri;
     for(int j=0;j<npri;j++){
        mo += moi[j]*coef[i_prim+j];
     }
    den += (nocc[i] * mo * mo);
  }
  return den;
}


double Field::DensitySYCL2(int norb, int npri,
                          int *icnt, int *vang,
                          double *r, double *coor, double *depris,
                          double *nocc, double *coef, double *moi){
  double den;
  double x = r[0];
  double y = r[1];
  double z = r[2];

  den=0.0;
  for(int i=0; i < norb; i++){
    double mo = 0.0;
     const int i_prim= i*npri;
    for(int j=0; j < npri; j++){
      const int vj = 3*j;
      const int centroj = 3*icnt[j];
      const double difx = x - coor[centroj];
      const double dify = y - coor[centroj+1];
      const double difz = z - coor[centroj+2];
      const double rr = difx*difx + dify*dify + difz*difz;

      const double expo = exp(-depris[j]*rr);
      const double lx   = vang[vj];
      const double ly   = vang[vj+1];
      const double lz   = vang[vj+2];
      const double facx = pow(difx,lx);
      const double facy = pow(dify,ly);
      const double facz = pow(difz,lz);

      mo += facx*facy*facz*expo * coef[i_prim+j];
    }
    den += nocc[i] * mo * mo;

  }

  return den;
}

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
  double delta = 0.5;

  int npoints_x = int((xmax - xmin) / delta);
  int npoints_y = int((ymax - ymin) / delta);
  int npoints_z = int((zmax - zmin) / delta);

  std::cout << " Points ( " << npoints_x << "," << npoints_y << "," << npoints_z
            << ")" << std::endl;
  std::cout << " TotalPoints : " << npoints_x * npoints_y * npoints_z
            << std::endl;

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
      cnt = 0;
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

void Field::evalDensity2D() {

  vector<double> field;
  double xmin = -10.0, xmax = 10.0;
  double ymin = -10.0, ymax = 10.0;
  double delta = 0.1;

  int npoints_x = int((xmax - xmin) / delta);
  int npoints_y = int((ymax - ymin) / delta);

  std::cout << " Points ( " << npoints_x << "," << npoints_y << ")  ";
  std::cout << " total points : " << npoints_x * npoints_y << std::endl;

  std::ofstream fout("density2d.dat");
  if (fout.is_open()) {

    for (int i = 0; i < npoints_x; i++) {
      double x = xmin + i * delta;
      for (int j = 0; j < npoints_y; j++) {
        double y = ymin + j * delta;

        Rvector r(x, y, 0.0);
        double den = Density(r);

        fout << std::setw(13) << std::setprecision(6) << std::fixed
             << x * 0.529177 << ' ';
        fout << std::setw(13) << std::setprecision(6) << std::fixed
             << y * 0.529177 << ' ';
        fout << std::setw(13) << std::setprecision(6) << std::fixed << den;
        fout << std::endl;
      }
      fout << std::endl;
    }
  }
}

void Field::evalDensity_sycl() {

  cl::sycl::queue q(cl::sycl::default_selector{});
  std::cout << " Running on " << q.get_device().get_info<cl::sycl::info::device::name>() << std::endl;

  double xmin = -10.0, xmax = 10.0;
  double ymin = -10.0, ymax = 10.0;
  double zmin = -5.0, zmax = 5.0;
  double delta = 0.5;
  vector<double> field;


  int npoints_x = int((xmax - xmin) / delta);
  int npoints_y = int((ymax - ymin) / delta);
  int npoints_z = int((zmax - zmin) / delta);
  const size_t nsize = npoints_x * npoints_y * npoints_z;
  int natm = wf.natm;
  int npri = wf.npri;
  int norb = wf.norb;
  int *icnt = new int[npri];
  int *vang = new int[3*npri];
  double *coor = new double [3*natm];
  double *eprim = new double [npri];
  double *moi = new double[npri];
  double *coef = new double [npri*norb];
  double *nocc = new double[norb];
  double *field_local = new double[nsize];

  std::cout << " Points ( " << npoints_x << "," << npoints_y << "," << npoints_z
            << ")" << std::endl;
  std::cout << " TotalPoints : " << npoints_x * npoints_y * npoints_z
            << std::endl;

  for(int i=0; i<npri;i++){
    icnt[i] = wf.icntrs[i];
    vang[3*i] = wf.vang[3*i];
    vang[3*i+1] = wf.vang[3*i+1];
    vang[3*i+2] = wf.vang[3*i+2];
    eprim[i] = wf.depris[i];
    moi[i] = 0.;
  }
  for(int i=0; i<norb; i++){
    nocc[i] =  wf.dnoccs[i];
    for(int j=0; j<npri; j++){
        coef[i*npri + j] = wf.dcoefs[i*npri + j];
    }
  }
  for(int i=0; i<natm; i++){
    Rvector R(wf.atoms[i].getCoors());
    coor[3*i] = R.get_x();
    coor[3*i+1] = R.get_y();
    coor[3*i+2] = R.get_z();
  }

  cl::sycl::buffer<int, 1>   icnt_buff   (icnt,  cl::sycl::range<1>(npri));
  cl::sycl::buffer<int, 1>   vang_buff   (vang,  cl::sycl::range<1>(3*npri));
  cl::sycl::buffer<double, 1> coor_buff  (coor,  cl::sycl::range<1>(3*natm));
  cl::sycl::buffer<double, 1> eprim_buff (eprim, cl::sycl::range<1>(npri));
  cl::sycl::buffer<double, 1> moi_buff   (moi,   cl::sycl::range<1>(npri));
  cl::sycl::buffer<double, 1> coef_buff  (coef,  cl::sycl::range<1>(npri*norb));
  cl::sycl::buffer<double, 1> nocc_buff  (nocc,  cl::sycl::range<1>(norb));
  cl::sycl::buffer<double, 1> field_buff (field_local, cl::sycl::range<1>(nsize));

  q.submit([&](cl::sycl::handler &h){
    auto field_acc = field_buff.get_access<cl::sycl::access::mode::write>(h);
    auto moi_acc = moi_buff.get_access<cl::sycl::access::mode::write>(h);
    auto icnt_acc = icnt_buff.get_access<cl::sycl::access::mode::read>(h);
    auto vang_acc = vang_buff.get_access<cl::sycl::access::mode::read>(h);
    auto coor_acc = coor_buff.get_access<cl::sycl::access::mode::read>(h);
    auto eprim_acc = eprim_buff.get_access<cl::sycl::access::mode::read>(h);
    auto coef_acc = coef_buff.get_access<cl::sycl::access::mode::read>(h);
    auto nocc_acc = nocc_buff.get_access<cl::sycl::access::mode::read>(h);

    h.parallel_for<class Field2>(cl::sycl::range<1>(nsize), [=](cl::sycl::id<1> idx){
      double cart[3];
      int k = (int) idx % npoints_z;
      int j = ((int) idx/npoints_z) % npoints_y;
      int i = (int) idx / (npoints_z * npoints_y);

      cart[0] = xmin + i * delta;
      cart[1] = ymin + j * delta;
      cart[2] = zmin + k * delta;
      int *icnt_ptr = icnt_acc.get_pointer();
      int *vang_ptr = vang_acc.get_pointer();
      double *coor_ptr = coor_acc.get_pointer();
      double *eprim_ptr = eprim_acc.get_pointer();
      double *nocc_ptr = nocc_acc.get_pointer();
      double *coef_ptr = coef_acc.get_pointer();
      double *moi_ptr = moi_acc.get_pointer();

      field_acc[idx] = DensitySYCL2(norb, npri, icnt_ptr, vang_ptr, cart, coor_ptr, eprim_ptr, nocc_ptr, coef_ptr, moi_ptr);
    });
  });
  q.wait();

  for(int i=0; i<nsize; i++)
     field.push_back(field_local[i]);

  dumpCube(xmin, ymin, zmin, delta, npoints_x, npoints_y, npoints_z, field);
  dumpXYZ();

  delete[] icnt;
  delete[] vang;
  delete[] coor;
  delete[] eprim;
  delete[] moi;
  delete[] coef;
  delete[] nocc;
  delete[] field_local;
}
