#include "propagate.hpp"

#define DEBUG

void propagate(Params * p) {
  // CVode variables
  void * cvode_mem = NULL;                      // pointer to block of CVode memory
  N_Vector y, yout;                     // arrays of populations

  int flag;
  realtype * dmt = NULL;                                // density matrix in time
  realtype * wfnt = NULL;                               // wave function in time
  realtype t0 = 0.0;                            // initial time
  realtype t = 0;
  realtype tret = 0;                                    // time returned by the solver


  if (p->wavefunction) {
    // Create the array to store the wavefunction in time
    wfnt = new realtype [2*p->NEQ*(p->numOutputSteps+1)];
    initializeArray(wfnt, 2*p->NEQ*(p->numOutputSteps+1), 0.0);
    for (int ii = 0; ii < 2*p->NEQ*(p->numOutputSteps+1); ii++) {
      wfnt[ii] = p->startWfn[ii];
    }
  }
  else {
    // Create the array to store the density matrix in time
    dmt = new realtype [2*p->NEQ2*(p->numOutputSteps+1)];
    initializeArray(dmt, 2*p->NEQ2*(p->numOutputSteps+1), 0.0);
    for (int ii = 0; ii < 2*p->NEQ2*(p->numOutputSteps+1); ii++) {
      dmt[ii] = p->startDM[ii];
    }
  }


  //// SET UP CVODE VARIABLES


#ifdef DEBUG
  std::cout << "\nCreating N_Vectors.\n";
  if (p->wavefunction) {
    std::cout << "\nProblem size is " << 2*p->NEQ << " elements.\n";
  }
  else {
    std::cout << "\nProblem size is " << 2*p->NEQ2 << " elements.\n";
  }
#endif
  // Creates N_Vector y with initial populations which will be used by CVode//
  if (p->wavefunction) {
    y = N_VMake_Serial(2*p->NEQ, &(p->startWfn[0]));
  }
  else {
    y = N_VMake_Serial(2*p->NEQ2, &(p->startDM[0]));
  }

  // the vector yout has the same dimensions as y
  yout = N_VClone(y);

  // Make outputs independent of time propagation ////////////////////////////

  computeGeneralOutputs(p);

  // create CVode object
  // this is a stiff problem, I guess?
#ifdef DEBUG
  std::cout << "\nCreating cvode_mem object.\n";
#endif
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  flag = CVodeSetUserData(cvode_mem, (void *) p);

#ifdef DEBUG
  std::cout << "\nInitializing CVode solver.\n";
#endif
  // initialize CVode solver //

  if (p->wavefunction) {
    //flag = CVodeInit(cvode_mem, &RHS_WFN, t0, y);
    flag = CVodeInit(cvode_mem, &RHS_WFN_SPARSE, t0, y);
  }
  else {
    if (p->kinetic) {
      flag = CVodeInit(cvode_mem, &RHS_DM_RELAX, t0, y);
    }
    else if (p->rta) {
      flag = CVodeInit(cvode_mem, &RHS_DM_RTA, t0, y);
      //flag = CVodeInit(cvode_mem, &RHS_DM_RTA_BLAS, t0, y);
    }
    else if (p->dephasing) {
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
  flag = CVodeSStolerances(cvode_mem, p->reltol, p->abstol);

#ifdef DEBUG
  std::cout << "\nAttaching linear solver module.\n";
#endif
  // attach linear solver module //
  if (p->wavefunction) {
    flag = CVDense(cvode_mem, 2*p->NEQ);
  }
  else {
    // Diagonal approximation to the Jacobian saves memory for large systems
    flag = CVDiag(cvode_mem);
  }

  //// CVODE TIME PROPAGATION


#ifdef DEBUG
  std::cout << "\nAdvancing the solution in time.\n";
#endif
  for (int ii = 1; ii <= p->numsteps; ii++) {
    t = (p->tout*((double) ii)/((double) p->numsteps));
    flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
    std::cout << std::endl << "CVode flag at step " << ii << ": " << flag << std::endl;
#endif
    if ((ii % (p->numsteps/p->numOutputSteps) == 0) || (ii == p->numsteps)) {
      // show progress in stdout
      if (p->progressStdout) {
        fprintf(stdout, "\r%-.2lf percent done", ((double)ii/((double)p->numsteps))*100);
        fflush(stdout);
      }
      // show progress in a file
      if (p->progressFile) {
        std::ofstream progressFile("progress.tmp");
        progressFile << ((double)ii/((double)p->numsteps))*100 << " percent done." << std::endl;
        progressFile.close();
      }
      if (p->wavefunction) {
        updateWfn(yout, wfnt, ii*p->numOutputSteps/p->numsteps, p);
      }
      else {
        updateDM(yout, dmt, ii*p->numOutputSteps/p->numsteps, p);
      }
    }
  }

  // Compute density outputs.
#ifdef DEBUG
  std::cout << "Computing outputs..." << std::endl;
#endif
  if (p->wavefunction) {
    computeWfnOutput(wfnt, p);
  }
  else {
    computeDMOutput(dmt, p);
  }
#ifdef DEBUG
  std::cout << "done computing outputs" << std::endl;
#endif

  // do analytical propagation
  if (p->analytical && (! p->bridge_on)) {
    computeAnalyticOutputs(p);
  }


  // clean up //////////////////////////////////////////////////////////////////


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
  if (p->wavefunction) {
    delete [] wfnt;
  }
  else {
    delete [] dmt;
  }

  return;
}