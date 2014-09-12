#include "rhs.hpp"

// #define DEBUG_RHS
//
// DEBUGf flag: general output at each CVode step
// #define DEBUGf
//
// DEBUGf_DM flag: DEBUGf for density matrix EOM
// #define DEBUGf_DM

// #define DEBUG_TORSION
// #define DEBUG_DYNAMIC_MU


/* Right-hand-side equation for wavefunction */
int RHS_WFN(realtype t, N_Vector y, N_Vector ydot, void * data) {
  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
#ifdef DEBUG
    std::cout << "Updating Hamiltonian at time " << t << std::endl;
#endif
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }
  else {
#ifdef DEBUG
    std::cout << "Not updating Hamiltonian at time " << t << std::endl;
#endif
  }

  // extract parameters from p
  realtype * H = &(p->H)[0];
  //realtype * H = &(p->H_lo)[0];
  int N = p->NEQ;

  // get pointer to y, ydot data
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

#ifdef DEBUG_RHS
  std::cout << "Pointer to ydot is " << ydotp << std::endl;
#endif

  // set up BLAS variables
  // const char TRANS = 'n';
  const char UPLO = 'l';
  double beta = 0.0;
  double alpha_re = 1.0;  // alpha value for real part of wfn derivative
  double alpha_im = -1.0; // alpha value for imag part of wfn derivative
  int inc = 1;

  // Re(\dot{\psi}) = \hat{H}Im(\psi)
  //DGEMV(&TRANS, &N, &N, &alpha_re, &H[0], &N, &yp[N], &inc, &beta, &ydotp[0], &inc);
  DSYMV(&UPLO, &N, &alpha_re, &H[0], &N, &yp[N], &inc, &beta, &ydotp[0], &inc);

  // Im(\dot{\psi}) = -i\hat{H}Re(\psi)
  //DGEMV(&TRANS, &N, &N, &alpha_im, &H[0], &N, &yp[0], &inc, &beta, &ydotp[N], &inc);
  DSYMV(&UPLO, &N, &alpha_im, &H[0], &N, &yp[0], &inc, &beta, &ydotp[N], &inc);

  return 0;
}

int RHS_WFN_SPARSE(realtype t, N_Vector y, N_Vector ydot, void * data) {
  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // extract parameters from p
  realtype * H = &(p->H_sp)[0];
  int * columns = &(p->H_cols)[0];
  int * rowind = &(p->H_rowind)[0];

  //realtype * H = &(p->H_lo)[0];
  int N = p->NEQ;

  // set up MKL variables
  char transa = 'n';
  double alpha_re = 1.0;  // alpha value for real part of wfn derivative
  double alpha_im = -1.0; // alpha value for imag part of wfn derivative
  double beta = 0.0;
  char matdescra [6] = {'S', // symmetric matrix
      'L', // lower triangle
      'N', // non-unit on diagonal
      'C', // zero-based indexing (C-style)
      '*', '*'}; // extra characters

  // Re(\dot{\psi}) = \hat{H}Im(\psi)
  mkl_dcsrmv(&transa, &N, &N, &alpha_re, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[N], &beta, &ydotp[0]);

  // Im(\dot{\psi}) = -i\hat{H}Re(\psi)
  mkl_dcsrmv(&transa, &N, &N, &alpha_im, &matdescra[0], &H[0], &columns[0],
             &rowind[0], &rowind[1], &yp[0], &beta, &ydotp[N]);

  return 0;
}

/* apply the kinetic relaxation model to one band of the system */
void RELAX_KINETIC(int bandFlag, realtype * yp, realtype * ydotp, Params * p) {
  double pop = 0.0;
  int start = bandStartIdx(bandFlag, p);
  int end = bandEndIdx(bandFlag, p);
  int Ni = bandNumStates(bandFlag, p);
  int N = p->NEQ;
  double g1;
  if (bandFlag == CONDUCTION) {
    g1 = p->gamma1;
  }
  else if (bandFlag == QD_CONDUCTION) {
    g1 = p->gamma1_c;
  }
  else {
    std::cerr << "WARNING [" << __FUNCTION__ << "]: unexpected bandFlag " << bandFlag << std::endl;
    std::cerr << "        setting g1 to 1.0" << std::endl;
    g1 = 1.0;
  }

  // don't do anything if there is no relaxation term
  if (g1 == 0.0) {
    return;
  }

  double mu = p->EF;
  double T = p->temperature;
  double * E = &(p->energies[start]);

  // sum current populations in band
  pop = 0.0;
  for (int ii = start; ii < end; ii++) {
    pop += yp[ii*N + ii];
  }

  // do the (N-1) pairs of states along diagonal
  double ePi, ePj;      // equilibrium populations, i and j indices
  double rel;           // relaxation term
  int Ii, Ij;           // indices

  // find equilibrium FDD
  double * fdd = new double [Ni];

  if ((p->dynamicMu) && (pop > 0.0)) {
    if (pop > 1.001) {  // test is against 1.001 since there may be some numerical drift
      std::cout << "WARNING [" << __FUNCTION__
        << "]: population in band is > 1; dynamic Fermi level may be spurious" << std::endl;
    }

    //// find bounds for Fermi level
    mu = findDynamicMu(pop, T, CONDUCTION, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "mu is " << mu << std::endl;
#endif
  }

  FDD(mu, T, fdd, E, Ni, pop);

  for (int ii = 0; ii < (Ni-1); ii++) {
    // precalculate indices and such
    Ii = (start + ii)*N + start + ii;
    Ij = Ii + N + 1;            // this index is the next diagonal element, so N+1 places up
    ePi = fdd[ii];
    ePj = fdd[ii+1];

    //// calculate contribution from relaxation

    // assuming \Gamma = k_{ij} + k_{ji}, "unimolecular" model
    rel = g1*(ePi*yp[Ij] - ePj*yp[Ii])/(ePi + ePj);

    // assuming \Gamma = k_{ij} + k_{ji}, "bimolecular" model
    // rel = g1*(yp[Ij]*(1-yp[Ii])*ePi*(1-ePj) - yp[Ii]*(1-yp[Ij])*ePj*(1-ePi))/(ePi+ePj-2*ePi*ePj);

    // assuming downward rates (j-->i) are the same, "bimolecular" model
    // rel = g1*(yp[Ij]*(1 - yp[Ii]) - ePj*(1-ePi)/(ePi*(1-ePj))*yp[Ii]*(1-yp[Ij]));

    // equal and opposite for the interaction of the two states
    ydotp[Ii] += rel;
    ydotp[Ij] -= rel;
  }

  delete [] fdd;

  return;
}

/* Right-hand-side equation for density matrix */
int RHS_DM_RELAX(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  double g2 = p->gamma2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // std::cout << "time is " << t << std::endl;
  // std::cout << ydotp << std::endl;
  // N_VPrint_Serial(y);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// relaxation in bulk conduction band

  if (p->kinetic) {
    RELAX_KINETIC(CONDUCTION, yp, ydotp, p);
  }

  //// relaxation in QD conduction band

  if (p->kineticQD) {
    RELAX_KINETIC(QD_CONDUCTION, yp, ydotp, p);
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
        //// real parts of ydot
        ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
        ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

        //// imaginary parts of ydot (lower triangle and complex conjugate)
        ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
        ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // dephasing
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUG_RHS
  std::cout << ydotp << " at time " << t << std::endl;
  N_VPrint_Serial(ydot);
#endif

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  dmf = fopen("dmf.out", "a");
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    fprintf(dmf, "\n");
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  fclose(dmf);
#endif

  return 0;
}

/* Right-hand-side equation for density matrix */
int RHS_DM_KINETIC(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  int Nk = p->Nk;
  int Ik = p->Ik;
  double * E = &(p->energies[Ik]);
  double g1 = p->gamma1;
  double g2 = p->gamma2;
  double T = p->temperature;
  double mu = p->EF;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// relaxation along diagonal
  // sum current populations in band
  double CBPop = 0.0;
  for (int ii = Ik; ii < (Ik + Nk); ii++) {
    CBPop += yp[ii*N + ii];
  }

  // do the (N-1) pairs of states along diagonal
  double ePi, ePj;      // equilibrium populations, i and j indices
  double rel;           // relaxation term
  int Ii, Ij;           // indices

  //// Kinetic relaxation model here

  // find equilibrium FDD
  double * fdd = new double [Nk];

  if ((p->dynamicMu) && (CBPop > 0.0)) {
    if (CBPop > 1.001) {        // test is against 1.001 since there may be some numerical drift
      std::cout << "WARNING [" << __FUNCTION__
        << "]: population in band is > 1; dynamic Fermi level may be spurious" << std::endl;
    }

    //// find bounds for Fermi level
    mu = findDynamicMu(CBPop, T, CONDUCTION, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "mu at time " << t << " is " << mu << std::endl;
#endif
  }

  FDD(mu, T, fdd, E, Nk, CBPop);

  for (int ii = 0; ii < (Nk-1); ii++) {
    // precalculate indices and such
    Ii = (Ik + ii)*N + Ik + ii;
    Ij = Ii + N + 1;            // this index is the next diagonal element, so N+1 places up
    ePi = fdd[ii];
    ePj = fdd[ii+1];

    // calculate contribution from relaxation
    rel = g1*(ePi*yp[Ij] - ePj*yp[Ii])/(ePi + ePj);

    // equal and opposite for the interaction of the two states
    ydotp[Ii] += rel;
    ydotp[Ij] -= rel;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
        //// real parts of ydot
        ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
        ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

        //// imaginary parts of ydot (lower triangle and complex conjugate)
        ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
        ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  dmf = fopen("dmf.out", "a");
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  fclose(dmf);
#endif

  // free fdd
  delete [] fdd;

  return 0;
}

/* Right-hand-side equation for density matrix */
int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "a");
#endif

  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

#ifdef DEBUG_RHS
  // print Hamiltonian
  std::cout << "Hamiltonian at time " << t << std::endl;
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      if (jj == 0) {
        fprintf(stdout, "%.9e", H[ii*N]);
      }
      else {
        fprintf(stdout, " %.9e", H[ii*N + jj]);
      }
    }
    fprintf(stdout, "\n");
  }
  std::cout << std::endl;
  // print DM
  std::cout << "DM at time " << t << std::endl;
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      if (jj == 0) {
        fprintf(stdout, "%.9e", yp[ii*N]);
      }
      else {
        fprintf(stdout, " %.9e", yp[ii*N + jj]);
      }
    }
    // Imaginary part
    fprintf(stdout, " ");
    for (int jj = 0; jj < N; jj++) {
      fprintf(stdout, " %.9e", yp[ii*N + jj + N2]);
    }
    fprintf(stdout, "\n");
  }
  std::cout << std::endl << std::endl;
#endif

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
        //// real parts of ydot
        ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
        ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

        //// imaginary parts of ydot (lower triangle and complex conjugate)
        ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
        ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  std::cout << "Closing output file for density matrix coefficients in time.\n";
  fclose(dmf);
#endif

#ifdef DEBUG_RHS
  // print DM'
  std::cout << "DM' at time " << t << " after propogator applied:" << std::endl;
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      if (jj == 0) {
        fprintf(stdout, "%.9e", ydotp[ii*N]);
      }
      else {
        fprintf(stdout, " %.9e", ydotp[ii*N + jj]);
      }
    }
    // Imaginary part
    fprintf(stdout, " ");
    for (int jj = 0; jj < N; jj++) {
      fprintf(stdout, " %.9e", ydotp[ii*N + jj + N2]);
    }
    fprintf(stdout, "\n");
  }
  std::cout << std::endl << std::endl;
#endif

// #ifdef DEBUG_RHS
//   std::cout << ydotp << " at time " << t << std::endl;
//   N_VPrint_Serial(ydot);
// #endif

  return 0;
}

/* Right-hand-side equation for density matrix using BLAS */
int RHS_DM_BLAS(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "a");
#endif

  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // extract parameters from p
  double * H = &(p->H)[0];
  int N = p->NEQ;
  int N2 = p->NEQ2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    // only update if at a new time point
    if ((t > 0.0) && (t != p->lastTime)) {
      updateHamiltonian(p, t);
      // update time point
      p->lastTime = t;
    }
  }

  // initialize ydot

#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  char LEFT = 'l';
  char RGHT = 'r';
  char UP = 'l';
  double ONE = 1.0;
  double NEG = -1.0;
  double ZERO = 0.0;

  // Re(\dot{\rho}) += H*Im(\rho)
  DSYMM(&LEFT, &UP, &N, &N, &ONE, &H[0], &N, &yp[N2], &N, &ZERO, &ydotp[0], &N);

  // Re(\dot{\rho}) -= Im(\rho)*H
  DSYMM(&RGHT, &UP, &N, &N, &NEG, &H[0], &N, &yp[N2], &N, &ONE, &ydotp[0], &N);

  // Im(\dot{\rho}) += i*Re(\rho)*H
  DSYMM(&RGHT, &UP, &N, &N, &ONE, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);

  // Im(\dot{\rho}) -= i*H*Re(\rho)
  DSYMM(&LEFT, &UP, &N, &N, &NEG, &H[0], &N, &yp[0], &N, &ONE, &ydotp[N2], &N);

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");

  std::cout << "Closing output file for density matrix coefficients in time.\n";
  fclose(dmf);
#endif

  return 0;
}

/* find Fermi level based on sum of population in band */
double findDynamicMu(double pop, double T, int bandFlag, Params * p) {
  double summ, lower, upper;

  // starting point
  double mu = 0.0;
  if (bandFlag == CONDUCTION) {
    mu = p->lastMu;
  }
  else if (bandFlag == QD_CONDUCTION) {
    mu = p->lastMuQD;
  }

  summ = FDDSum(mu, T, bandFlag, p);

  // just in case mu is exactly zero
  if (bool(fabs((summ - pop) > 1e-10))) {
    // otherwise search in increments of 1 a.u. energy
    if (summ < pop) {
      lower = -1.0;
      upper = 0.0;
      while (summ < pop) {
        lower++;
        upper++;
#ifdef DEBUG_DYNAMIC_MU
        std::cout << "Lower bound for mu: " << lower << std::endl;
        std::cout << "Upper bound for mu: " << upper << std::endl;
#endif
        summ = FDDSum(upper, T, bandFlag, p);
      }
    }
    else {
      lower = 0.0;
      upper = 1.0;
      while (summ > pop) {
        lower--;
        upper--;
#ifdef DEBUG_DYNAMIC_MU
        std::cout << "Lower bound for mu: " << lower << std::endl;
        std::cout << "Upper bound for mu: " << upper << std::endl;
#endif
        summ = FDDSum(lower, T, bandFlag, p);
      }
    }

    // do a binary search for mu
    mu = FDDBinarySearch(lower, upper, T, pop, bandFlag, p);
  }

  // store the mu for the next time step
  if (bandFlag == CONDUCTION) {
    p->lastMu = mu;
  }
  else if (bandFlag == QD_CONDUCTION) {
    p->lastMuQD = mu;
  }

  return mu;
}

/* Do a binary search to find the value of the Fermi level which makes the
 * sum of populations in a band add up to a certain value.
 */
double FDDBinarySearch(double lower, double upper, double T, double n,
    int bandFlag, Params * p) {
  double mid, summ;

  if (fabs(upper - lower) < 1e-10) {
    return lower;
  }
  else {
    mid = (upper + lower)/2.0;
    summ = FDDSum(mid, T, bandFlag, p);
#ifdef DEBUG_DYNAMIC_MU
    std::cout << "Binary search lower bound: " << lower << std::endl;
    std::cout << "Binary search upper bound: " << upper << std::endl;
    std::cout << "Binary search middle: " << mid << std::endl;
    std::cout << "Binary search target value: " << n << std::endl;
    std::cout << "Binary search summ value: " << summ << std::endl;
#endif
    if (summ > n) {
#ifdef DEBUG_DYNAMIC_MU
      std::cout << "value is in lower half of bounds" << std::endl;
      std::cout << std::endl;
#endif
      return FDDBinarySearch(lower, mid, T, n, bandFlag, p);
    }
    else {
#ifdef DEBUG_DYNAMIC_MU
      std::cout << "value is in upper half of bounds" << std::endl;
      std::cout << std::endl;
#endif
      return FDDBinarySearch(mid, upper, T, n, bandFlag, p);
    }
  }
}

/* Add up the populations in a band with a Fermi-Dirac distribution of population
 */
double FDDSum(double mu, double T, int bandFlag, Params * p) {
  double summ = 0.0;
  int start = bandStartIdx(bandFlag, p);
  int end = bandEndIdx(bandFlag, p);
  double beta = 3.185e5/T;
  double * E = &(p->energies[start]);

  for (int ii = start; ii < end; ii++) {
    summ += 1.0/(1.0 + exp((E[ii] - mu)*beta));
  }

  return summ;
}

/* fills the array fdd with Fermi-Dirac populations, normalized to a population
 * P.
 */
void FDD(double mu, double T, double * fdd, double * E, int N, double P) {
  double beta = 3.185e5/T;

  for (int ii = 0; ii < N; ii++) {
    fdd[ii] = 1.0/(1.0 + exp((E[ii] - mu)*beta));
  }

  if (P != 0.0) {
    double norm = P/std::accumulate(fdd, fdd+N, 0.0);

    for (int ii = 0; ii < N; ii++) {
      fdd[ii] *= norm;
    }
  }

  return;
}

/* Right-hand-side equation for density matrix
 * using dephasing */
int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "a");
#endif


  // data is a pointer to the params struct
  Params * p;
  p = (Params *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g2 = p->gamma2;

  // more compact notation for N_Vectors
  realtype * yp = N_VGetArrayPointer(y);
  realtype * ydotp = N_VGetArrayPointer(ydot);

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    updateHamiltonian(p, t);
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    ydotp[ii] = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      ydotp[ii*N + ii] += 2*H[ii*N + jj]*yp[jj*N + ii + N2];
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
        //// real parts of ydot
        ydotp[ii*N + jj] += H[ii*N + kk]*yp[kk*N + jj + N2];
        ydotp[ii*N + jj] -= yp[ii*N + kk + N2]*H[kk*N + jj];

        //// imaginary parts of ydot (lower triangle and complex conjugate)
        ydotp[ii*N + jj + N2] -= H[ii*N + kk]*yp[kk*N + jj];
        ydotp[ii*N + jj + N2] += yp[ii*N + kk]*H[kk*N + jj];
      }
      // the complex conjugate
      ydotp[jj*N + ii] = ydotp[ii*N + jj];
      ydotp[jj*N + ii + N2] = -1*ydotp[ii*N + jj + N2];

      // relaxation
      ydotp[ii*N + jj] -= g2*yp[ii*N + jj];
      ydotp[ii*N + jj + N2] -= g2*yp[ii*N + jj + N2];
      ydotp[jj*N + ii] -= g2*yp[jj*N + ii];
      ydotp[jj*N + ii + N2] -= g2*yp[jj*N + ii + N2];
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", ydotp[ii*N + jj], ydotp[ii*N + jj + N2]);
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}