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

/* gives the equilibrated FDD for the system */
void buildFDD(struct PARAMETERS * p, N_Vector y, std::vector<double> & fdd) {
  //// "fine structure constant" -- conversion from index to wave vector
  double x2 = 8e-8;

  //// calculate n_e and e_kin
  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(x2,3));
  for (int ii = 0; ii < p->Nk; ii++) {
    ne += factor*pow(ii,2)*NV_Ith_S(y, ii*p->NEQ + ii);
    ekin += factor*pow(ii,4)*NV_Ith_S(y, ii*p->NEQ + ii)/(2*p->me*pow(x2,2));
  }
  
  //// find the inverse temperature (beta)
  int iter = 0;
  int maxiter = 200;
  double tol = 1e-12;
  double K1 = 4.8966851;		// constants
  double K2 = 0.04496457;
  double K3 = 0.133376;
  double X = 4*p->me*pow(M_PI/(2*p->me),1.5);
  double bn = 1.0;		// intermediate values of beta; bn is higher iteration
  double bm = 0.0;
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
  // loop applies Newton-Rhapson method to get zero of function
  while ((fabs(bn - bm) > tol) && (iter < maxiter)) {
    bm = bn;
    f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
    fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*pow(bm,0.5));
    bn = bm - f/fp;
    iter++;
    std::cout << "Iteration " << iter << std::endl;
    std::cout << "f  " << f << std::endl;
    std::cout << "f' " << fp << std::endl;
  }
  std::cout << std::endl;

  //// use beta to find chemical potential
  double mue = 0.0;
  mue = pow(3*ne,2.0/3.0)*pow(M_PI,4.0/3.0)/(2*p->me) - 2*pow(M_PI,2.0/3.0)*p->me/(pow(3*ne,2.0/3.0)*6*pow(bn,2));

  // TODO account for temperature dropping in time
  std::cout << "inverse temp is " << bn << std::endl;
  std::cout << std::endl;
  for (int ii = 0; ii < p->Nk; ii++) {
    fdd[ii] = 1.0/(1.0 + exp((ekin - mue)*bn));
    std::cout << 1.0/(1.0 + exp((ekin - mue)*bn)) << " ";
    std::cout << "FDD[" << ii << "]: " << fdd[ii] << std::endl;
  }

  return;
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
  realtype g1 = p->gamma1;
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
  //   get equilibrium FDD populations
  std::vector<double> fdd(p->Nk);
  buildFDD(p, y, fdd);

#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      NV_Ith_S(ydot, ii*N + ii) += 2*H[ii*N + jj]*NV_Ith_S(y, jj*N + ii + N2)
	- g1*(NV_Ith_S(y, ii*N + ii) - fdd[ii]);
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

