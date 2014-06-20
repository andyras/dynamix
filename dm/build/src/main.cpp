#include <getopt.h>

#include "dynamix.hpp"
#include "propagate.hpp"

int main (int argc, char * argv[]) {
  // Declare variables

  // Struct of parameters
  Params p;
  time_t startRun;                              // time at start of log
  time_t endRun;                                        // time at end of log
  struct tm * currentTime = NULL;                       // time structure for localtime
  FILE * log;                                   // log file with run times

  double summ = 0;                      // sum variable

  // ---- process command line flags ---- //
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
          std::cerr << "ERROR: option -i requires argument ("
                    << insDir << ") to have a trailing slash (/)." << std::endl;
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

  // print start time to log file //////////////////////////////////////////////

  if (isOutput(p.outs, "log.out")) {
    log = fopen("log.out", "w"); // note that this file is closed at the end of the program
  }
  time(&startRun);
  currentTime = localtime(&startRun);
  if (isOutput(p.outs, "log.out")) {
    fprintf(log, "Run started at %s\n", asctime(currentTime));
  }

  if (isOutput(p.outs, "log.out") && p.laser_on) {
    // make a note about the laser intensity.
    fprintf(log,"The laser intensity is %.5e W/cm^2.\n\n",pow(p.pumpAmpl,2)*3.5094452e16);
  }

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  if (isOutput(p.outs, "log.out")) {
    // make a note in the log about system timescales
    double tau = 0;             // fundamental system timescale
    if (p.Nk == 1) {
      fprintf(log, "\nThe timescale (tau) is undefined (Nk == 1).\n");
    }
    else {
      if (p.bridge_on) {
        if (p.scale_bubr) {
          tau = 1.0/(2*p.Vbridge[0]*M_PI);
        }
        else {
          tau = ((p.kBandTop - p.kBandEdge)/(p.Nk - 1))/(2*pow(p.Vbridge[0],2)*M_PI);
        }
      }
      else {
        if (p.scale_buqd) {
          tau = 1.0/(2*p.Vnobridge[0]*M_PI);
        }
        else {
          tau = ((p.kBandTop - p.kBandEdge)/(p.Nk - 1))/(2*pow(p.Vnobridge[0],2)*M_PI);
        }
      }
      fprintf(log, "\nThe timescale (tau) is %.9e a.u.\n", tau);
    }
  }

  // Make plot files ///////////////////////////////////////////////////////////
  makePlots(&p);

  // only do propagation if not just making plots //////////////////////////////

  if (! p.justPlots) {
    propagate(&p);

    // finalize log file //
    time(&endRun);
    currentTime = localtime(&endRun);
    if (isOutput(p.outs, "log.out")) {
      fprintf(log, "Run ended at %s\n", asctime(currentTime));
      fprintf(log, "Run took %.3g seconds.\n", difftime(endRun, startRun));
      fclose(log);                                      // note that the log file is opened after variable declaration
    }
    if (p.progressStdout) {
      printf("\nRun took %.3g seconds.\n", difftime(endRun, startRun));
    }
  }

  std::cout << "whoo" << std::endl;

  return 0;
}