#include "rhs.hpp"

//#define DEBUG_RHS

/* Updates the Hamiltonian with the time-dependent torsional coupling */
void updateTorsionV(PARAMETERS * p, realtype t) {
  double torsionValue = p->torsionV->value(t);

  // bridge is off, coupling is between k and c states
  if (!(p->bridge_on)) {
#ifdef DEBUG_RHS
    std::cout << "torsion between k and c states" << std::endl;
#endif
    for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
      for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
	p->H[ii*p->NEQ + jj] = torsionValue;
	p->H[jj*p->NEQ + ii] = torsionValue;
      }
    }
  }
  // torsion is at first bridge coupling
  else if (p->torsionSite == 0) {
#ifdef DEBUG_RHS
    std::cout << "torsion between k states and bridge" << std::endl;
#endif
    for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
      p->H[ii*p->NEQ + p->Ib] = torsionValue;
      p->H[p->Ib*p->NEQ + ii] = torsionValue;
    }
  }
  // torsion is at last bridge coupling
  else if (p->torsionSite == p->Nb) {
#ifdef DEBUG_RHS
    std::cout << "torsion between bridge and c states" << std::endl;
#endif
    for (int ii = p->Ic; ii < (p->Ic + p->Nc); ii++) {
      p->H[ii*p->NEQ + p->Ib + p->Nb - 1] = torsionValue;
      p->H[(p->Ib + p->Nb - 1)*p->NEQ + ii] = torsionValue;
    }
  }
  // torsion is between bridge sites
  else {
#ifdef DEBUG_RHS
    std::cout << "torsion between bridge sites " << p->torsionSite - 1
              << " and " << p->torsionSite << "." << std::endl;
#endif
    p->H[(p->Ib + p->torsionSite - 1)*p->NEQ + p->Ib + p->torsionSite] = torsionValue;
    p->H[(p->Ib + p->torsionSite)*p->NEQ + p->Ib + p->torsionSite - 1] = torsionValue;
  }

  return;
}

/* Right-hand-side equation for density matrix */
int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  
  // if torsion is on, update Hamiltonian
  if (p->torsion) {
    updateTorsionV(p, t);
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    NV_Ith_S(ydot, ii) = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      NV_Ith_S(ydot, ii*N + ii) += 2*H[ii*N + jj]*NV_Ith_S(y, jj*N + ii + N2);
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	NV_Ith_S(ydot, ii*N + jj) += H[ii*N + kk]*NV_Ith_S(y, kk*N + jj + N2);
	NV_Ith_S(ydot, ii*N + jj) -= NV_Ith_S(y, ii*N + kk + N2)*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	NV_Ith_S(ydot, ii*N + jj + N2) -= H[ii*N + kk]*NV_Ith_S(y, kk*N + jj);
	NV_Ith_S(ydot, ii*N + jj + N2) += NV_Ith_S(y, ii*N + kk)*H[kk*N + jj];
      }
      // the complex conjugate
      NV_Ith_S(ydot, jj*N + ii) = NV_Ith_S(ydot, ii*N + jj);
      NV_Ith_S(ydot, jj*N + ii + N2) = -1*NV_Ith_S(ydot, ii*N + jj + N2);
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", NV_Ith_S(ydot, ii*N + jj), NV_Ith_S(ydot, ii*N + jj + N2));
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}

/* Right-hand-side equation for density matrix
 * using relaxation time approximation (RTA) */
int RHS_DM_RTA(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g2 = p->gamma2;

  // if torsion is on, update Hamiltonian
  if (p->torsion) {
    updateTorsionV(p, t);
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    NV_Ith_S(ydot, ii) = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      NV_Ith_S(ydot, ii*N + ii) += 2*H[ii*N + jj]*NV_Ith_S(y, jj*N + ii + N2);
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	NV_Ith_S(ydot, ii*N + jj) += H[ii*N + kk]*NV_Ith_S(y, kk*N + jj + N2);
	NV_Ith_S(ydot, ii*N + jj) -= NV_Ith_S(y, ii*N + kk + N2)*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	NV_Ith_S(ydot, ii*N + jj + N2) -= H[ii*N + kk]*NV_Ith_S(y, kk*N + jj);
	NV_Ith_S(ydot, ii*N + jj + N2) += NV_Ith_S(y, ii*N + kk)*H[kk*N + jj];
      }
      // the complex conjugate
      NV_Ith_S(ydot, jj*N + ii) = NV_Ith_S(ydot, ii*N + jj);
      NV_Ith_S(ydot, jj*N + ii + N2) = -1*NV_Ith_S(ydot, ii*N + jj + N2);

      // relaxation
      NV_Ith_S(ydot, ii*N + jj) -= g2*NV_Ith_S(y, ii*N + jj);
      NV_Ith_S(ydot, ii*N + jj + N2) -= g2*NV_Ith_S(y, ii*N + jj + N2);
      NV_Ith_S(ydot, jj*N + ii) -= g2*NV_Ith_S(y, jj*N + ii);
      NV_Ith_S(ydot, jj*N + ii + N2) -= g2*NV_Ith_S(y, jj*N + ii + N2);
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", NV_Ith_S(ydot, ii*N + jj), NV_Ith_S(ydot, ii*N + jj + N2));
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}

/* Right-hand-side equation for density matrix
 * using dephasing */
int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data) {

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g2 = p->gamma2;

  // if torsion is on, update Hamiltonian
  if (p->torsion) {
    updateTorsionV(p, t);
  }

  // initialize ydot
  // THIS NEEDS TO BE HERE FOR SOME REASON EVEN IF ALL ELEMENTS ARE ASSIGNED LATER
#pragma omp parallel for
  for (int ii = 0; ii < 2*N2; ii++) {
    NV_Ith_S(ydot, ii) = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      NV_Ith_S(ydot, ii*N + ii) += 2*H[ii*N + jj]*NV_Ith_S(y, jj*N + ii + N2);
    }
  }

  //// off-diagonal
#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < ii; jj++) {
      for (int kk = 0; kk < N; kk++) {
	//// real parts of ydot
	NV_Ith_S(ydot, ii*N + jj) += H[ii*N + kk]*NV_Ith_S(y, kk*N + jj + N2);
	NV_Ith_S(ydot, ii*N + jj) -= NV_Ith_S(y, ii*N + kk + N2)*H[kk*N + jj];

	//// imaginary parts of ydot (lower triangle and complex conjugate)
	NV_Ith_S(ydot, ii*N + jj + N2) -= H[ii*N + kk]*NV_Ith_S(y, kk*N + jj);
	NV_Ith_S(ydot, ii*N + jj + N2) += NV_Ith_S(y, ii*N + kk)*H[kk*N + jj];
      }
      // the complex conjugate
      NV_Ith_S(ydot, jj*N + ii) = NV_Ith_S(ydot, ii*N + jj);
      NV_Ith_S(ydot, jj*N + ii + N2) = -1*NV_Ith_S(ydot, ii*N + jj + N2);

      // relaxation
      NV_Ith_S(ydot, ii*N + jj) -= g2*NV_Ith_S(y, ii*N + jj);
      NV_Ith_S(ydot, ii*N + jj + N2) -= g2*NV_Ith_S(y, ii*N + jj + N2);
      NV_Ith_S(ydot, jj*N + ii) -= g2*NV_Ith_S(y, jj*N + ii);
      NV_Ith_S(ydot, jj*N + ii + N2) -= g2*NV_Ith_S(y, jj*N + ii + N2);
    }
  }

#ifdef DEBUGf_DM
  fprintf(dmf, "%+.7e", t);
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      fprintf(dmf, " (%+.2e,%+.2e)", NV_Ith_S(ydot, ii*N + jj), NV_Ith_S(ydot, ii*N + jj + N2));
    }
  }
  fprintf(dmf, "\n");
#endif

  return 0;
}

