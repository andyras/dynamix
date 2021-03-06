#include <getopt.h>

#include "dynamix.hpp"
#include "log.hpp"
#include "propagate.hpp"

  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;

int main (int argc, char * argv[]) {
  // Set up logging ////////////////////////////////////////////////////////////

  initLog();

  // Struct of parameters //////////////////////////////////////////////////////

  Params p;

  // process command line flags ////////////////////////////////////////////////

  opterr = 0;
  int c;
  std::string insDir;
  /* process command line options */
  while ((c = getopt(argc, argv, "i:o:")) != -1) {
    switch (c) {
      case 'i':
        // check that it ends in a slash
        insDir = optarg;
        if (strcmp(&(insDir.at(insDir.length() - 1)), "/")) {
          BOOST_LOG_SEV(lg, error) << "option -i requires argument ("
                    << insDir << ") to have a trailing slash (/).";
          return 1;
        }
        else {
          // ---- assign input files ---- //
          p.inputFile = insDir + "parameters.in";
          p.cEnergiesInput = insDir + "c_energies.in";
          p.bEnergiesInput = insDir + "b_energies.in";
          p.VNoBridgeInput = insDir + "Vnobridge.in";
          p.VBridgeInput = insDir + "Vbridge.in";
        }
        break;
      case 'o':
        p.outputDir = optarg;
        break;
      case '?':
        if (optopt == 'i') {
          fprintf(stderr, "Option -%c requires a directory argument.\n", optopt);
        }
        else if (isprint(optopt)) {
          fprintf(stderr, "Unknown option -%c.\n", optopt);
        }
        else {
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        }
        return 1;
      default:
        continue;
    }
  }

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p.inputFile.c_str(), &p);

  initialize(&p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
#ifdef __USE_MKL__
  mkl_set_num_threads(p.nproc);
#endif

  // Make plot files ///////////////////////////////////////////////////////////

  makePlots(&p);

  // only do propagation if not just making plots //////////////////////////////

  if (! p.justPlots) {
    propagate(&p);
  }

  BOOST_LOG_SEV(lg, info) << "whoo";

  return 0;
}