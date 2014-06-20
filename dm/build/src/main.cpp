#include "dynamix.hpp"
#include <getopt.h>

int main (int argc, char * argv[]) {
  // Declare variables

  // Struct of parameters
  Params p;
  // CVode variables
  void * cvode_mem = NULL;                      // pointer to block of CVode memory
  N_Vector y, yout;                     // arrays of populations

  int flag;
  realtype * ydata = NULL;                              // pointer to ydata (contains all populations)
  realtype * dmt = NULL;                                // density matrix in time
  realtype * wfnt = NULL;                               // wave function in time
  realtype t0 = 0.0;                            // initial time
  realtype t = 0;
  realtype tret = 0;                                    // time returned by the solver
  time_t startRun;                              // time at start of log
  time_t endRun;                                        // time at end of log
  struct tm * currentTime = NULL;                       // time structure for localtime
#ifdef DEBUG
  FILE * realImaginary;                         // file containing real and imaginary parts of the wavefunction
#endif
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

  if (fileExists(p.outputDir)) {
    flag = mkdir(p.outputDir.c_str(), 0755);
  }
  // ---- TODO create output directory if it does not exist ---- //

  // assign parameters from input file /////////////////////////////////////////

  assignParams(p.inputFile.c_str(), &p);

  // Decide which output files to make /////////////////////////////////////////

#ifdef DEBUG
  std::cout << "Assigning outputs as specified in " << p.inputFile << "\n";
#endif
  assignOutputs(p.inputFile.c_str(), p.outs, &p);

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

  // TODO replace with initialize function
  initHamiltonian(&p);
  initWavefunction(&p);

  // set number of processors for OpenMP ///////////////////////////////////////

  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);

  // Assign coupling matrix ////////////////////////////////////////////////////

  p.buildCoupling();

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

  if (p.wavefunction) {
    // Create the array to store the wavefunction in time
    wfnt = new realtype [2*p.NEQ*(p.numOutputSteps+1)];
    initializeArray(wfnt, 2*p.NEQ*(p.numOutputSteps+1), 0.0);
  }
  else {
    // Create the array to store the density matrix in time
    dmt = new realtype [2*p.NEQ2*(p.numOutputSteps+1)];
    initializeArray(dmt, 2*p.NEQ2*(p.numOutputSteps+1), 0.0);
  }

#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif
  p.buildHamiltonian();


  //// SET UP CVODE VARIABLES


#ifdef DEBUG
  std::cout << "\nCreating N_Vectors.\n";
  if (p.wavefunction) {
    std::cout << "\nProblem size is " << 2*p.NEQ << " elements.\n";
  }
  else {
    std::cout << "\nProblem size is " << 2*p.NEQ2 << " elements.\n";
  }
#endif
  // Creates N_Vector y with initial populations which will be used by CVode//
  if (p.wavefunction) {
    y = N_VMake_Serial(2*p.NEQ, &(p.startWfn[0]));
  }
  else {
    y = N_VMake_Serial(2*p.NEQ2, &(p.startDM[0]));
  }

  // put in t = 0 information
  if (! p.wavefunction) {
    updateDM(y, dmt, 0, &p);
  }
  else {
    updateWfn(y, wfnt, 0, &p);
  }
  // the vector yout has the same dimensions as y
  yout = N_VClone(y);

#ifdef DEBUG
  realImaginary = fopen("real_imaginary.out", "w");
#endif

  // Make plot files
  makePlots(&p);

  // only do propagation if not just making plots
  if (! p.justPlots) {
    // Make outputs independent of time propagation
    computeGeneralOutputs(&p);
std::cout << "\n\n\nWHOOOOOT\n\n\n";

    // create CVode object
    // this is a stiff problem, I guess?
#ifdef DEBUG
    std::cout << "\nCreating cvode_mem object.\n";
#endif
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    flag = CVodeSetUserData(cvode_mem, (void *) &p);

#ifdef DEBUG
    std::cout << "\nInitializing CVode solver.\n";
#endif
    // initialize CVode solver //

    if (p.wavefunction) {
      //flag = CVodeInit(cvode_mem, &RHS_WFN, t0, y);
      flag = CVodeInit(cvode_mem, &RHS_WFN_SPARSE, t0, y);
    }
    else {
      if (p.kinetic) {
        flag = CVodeInit(cvode_mem, &RHS_DM_RELAX, t0, y);
      }
      else if (p.rta) {
        flag = CVodeInit(cvode_mem, &RHS_DM_RTA, t0, y);
        //flag = CVodeInit(cvode_mem, &RHS_DM_RTA_BLAS, t0, y);
      }
      else if (p.dephasing) {
        flag = CVodeInit(cvode_mem, &RHS_DM_dephasing, t0, y);
      }
      else {
        //flag = CVodeInit(cvode_mem, &RHS_DM, t0, y);
        flag = CVodeInit(cvode_mem, &RHS_DM_BLAS, t0, y);
      }
    }

#ifdef DEBUG
    std::cout << "\nSpecifying integration tolerances.\n";
#endif
    // specify integration tolerances //
    flag = CVodeSStolerances(cvode_mem, p.reltol, p.abstol);

#ifdef DEBUG
    std::cout << "\nAttaching linear solver module.\n";
#endif
    // attach linear solver module //
    if (p.wavefunction) {
      flag = CVDense(cvode_mem, 2*p.NEQ);
    }
    else {
      // Diagonal approximation to the Jacobian saves memory for large systems
      flag = CVDiag(cvode_mem);
    }

    //// CVODE TIME PROPAGATION


#ifdef DEBUG
    std::cout << "\nAdvancing the solution in time.\n";
#endif
    for (int ii = 1; ii <= p.numsteps; ii++) {
      t = (p.tout*((double) ii)/((double) p.numsteps));
      flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
      std::cout << std::endl << "CVode flag at step " << ii << ": " << flag << std::endl;
#endif
      if ((ii % (p.numsteps/p.numOutputSteps) == 0) || (ii == p.numsteps)) {
        // show progress in stdout
        if (p.progressStdout) {
          fprintf(stdout, "\r%-.2lf percent done", ((double)ii/((double)p.numsteps))*100);
          fflush(stdout);
        }
        // show progress in a file
        if (p.progressFile) {
          std::ofstream progressFile("progress.tmp");
          progressFile << ((double)ii/((double)p.numsteps))*100 << " percent done." << std::endl;
          progressFile.close();
        }
        if (p.wavefunction) {
          updateWfn(yout, wfnt, ii*p.numOutputSteps/p.numsteps, &p);
        }
        else {
          updateDM(yout, dmt, ii*p.numOutputSteps/p.numsteps, &p);
        }
      }
    }

#ifdef DEBUG
    fclose(realImaginary);
#endif

    // finalize log file //
    time(&endRun);
    currentTime = localtime(&endRun);
    if (isOutput(p.outs, "log.out")) {
      fprintf(log, "Final status of 'flag' variable: %d\n\n", flag);
      fprintf(log, "Run ended at %s\n", asctime(currentTime));
      fprintf(log, "Run took %.3g seconds.\n", difftime(endRun, startRun));
      fclose(log);                                      // note that the log file is opened after variable declaration
    }
    if (p.progressStdout) {
      printf("\nRun took %.3g seconds.\n", difftime(endRun, startRun));
    }

    // Compute density outputs.
#ifdef DEBUG
    std::cout << "Computing outputs..." << std::endl;
#endif
    if (p.wavefunction) {
      computeWfnOutput(wfnt, &p);
    }
    else {
      computeDMOutput(dmt, &p);
    }
#ifdef DEBUG
    std::cout << "done computing outputs" << std::endl;
#endif

    // do analytical propagation
    if (p.analytical && (! p.bridge_on)) {
      computeAnalyticOutputs(&p);
    }
  }


  //// CLEAN UP


#ifdef DEBUG
  fprintf(stdout, "Deallocating N_Vectors.\n");
#endif
  // deallocate memory for N_Vectors //
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(yout);

#ifdef DEBUG
  fprintf(stdout, "Freeing CVode memory.\n");
#endif
  // free solver memory //
  CVodeFree(&cvode_mem);

#ifdef DEBUG
  fprintf(stdout, "Freeing memory in main.\n");
#endif
  // delete all these guys
  if (p.wavefunction) {
    delete [] wfnt;
  }
  else {
    delete [] dmt;
  }

  std::cout << "whoo" << std::endl;

  return 0;
}