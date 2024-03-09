#include "Field.hpp"
#include "Wavefunction.hpp"
#include "version.hpp"

int main(int argc, char *argv[]) {
  std::cout << "Version: " << PROJECT_VER << std::endl;
  std::cout << "Compilation Date: " << __DATE__ << "  " << __TIME__
            << std::endl;
  std::cout << "Git SHA1: " << GIT_SHA1 << std::endl;

  Wavefunction wf;
  wf.loadWF(argv[1]);

  Field field(wf);

  field.evalDensity();

  // wf.printWF();
}
