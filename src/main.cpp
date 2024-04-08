#include "Field.hpp"
#include "Wavefunction.hpp"
#include "version.hpp"
#include <cstdlib>

int main(int argc, char *argv[]) {
  std::cout << "Version: " << PROJECT_VER << std::endl;
  std::cout << "Compilation Date: " << __DATE__ << "  " << __TIME__
            << std::endl;
  std::cout << "Git SHA1: " << GIT_SHA1 << std::endl;

  Wavefunction wf;
  if( argc != 2){
    std::cout << " We need more arguments try with:" << std::endl;
    std::cout << " ./" << argv[0] << " foo.wfx" << std::endl;
    exit(EXIT_FAILURE);
  }
  wf.loadWF(argv[1]);

  Field field(wf);

  field.evalDensity_sycl();

  // wf.printWF();
  exit(EXIT_SUCCESS);
}
