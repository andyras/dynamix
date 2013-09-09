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
  std::cout << "p->X2   " << p->X2 << std::endl;

  //// calculate n_e and e_kin
  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(p->X2,3));
  std::cout << "factor   " << factor << std::endl;
  for (int ii = 0; ii < p->Nk; ii++) {
    ne += factor*pow(ii,2)*NV_Ith_S(y, ii*p->NEQ + ii);
    std::cout << "population " << NV_Ith_S(y, ii*p->NEQ + ii) << std::endl;
    ekin += factor*pow(ii,4)*NV_Ith_S(y, ii*p->NEQ + ii)/(2*p->me*pow(p->X2,2));
  }
  std::cout << "ne   " << ne << std::endl;
  std::cout << "ekin " << ekin << std::endl;
  
  //// find the inverse temperature (beta)
  int iter = 0;
  int maxiter = 200;
  double tol = 1e-18;
  double K1 = 4.8966851;		// constants
  double K2 = 0.04496457;
  double K3 = 0.133376;
  double X = 4*ne*pow(M_PI/(2*p->me),1.5);
  double bn = 1e-10;		// intermediate values of beta; bn is higher iteration
  double bm = 0e-10;
  /*
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
  // loop applies Newton-Raphson method to get zero of function
  std::cout << "Newton-Raphson to find inverse temperature" << std::endl;
  std::cout << "X " << X << std::endl;
  while ((fabs(bn - bm) > tol) && (iter < maxiter)) {
    bm = bn;
    f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
    fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*pow(bm,0.5));
    bn = bm - f/fp;
    iter++;
    std::cout << "Iteration " << iter << std::endl;
    std::cout << "f  " << f << std::endl;
    std::cout << "f' " << fp << std::endl;
    std::cout << "bn " << bn << std::endl;
    std::cout << "bm " << bm << std::endl;
  }
  std::cout << std::endl;
  */
  
  double high = 1e10;
  double low = 1.0e-100;	// zero is a no-no because the function is a log
  double newVal = 0.5;		// intermediate value

  // check that f(high) and f(low) have opposite sign
  if (sgn<double>(b13(low, ekin, ne, K1, K2, K3, X))*sgn<double>(b13(high, ekin, ne, K1, K2, K3, X) > 0)) {
    std::cout << "ERROR: f(high) and f(low) have same sign!!!!" << std::endl;
    std::cout << "f(" << high << "): " << b13(high, ekin, ne, K1, K2, K3, X) << std::endl;
    std::cout << "f(" << low << "): " << b13(low, ekin, ne, K1, K2, K3, X) << std::endl;
  }

  // loop does binary search to find zero of function
  while ((iter < maxiter) && ((high-low) > tol)) {
    std::cout << "---ITERATION " << iter << "---" << std::endl;
    newVal = (high + low)/2.0;
    if (sgn<double>(b13(newVal, ekin, ne, K1, K2, K3, X)) == sgn<double>(b13(low, ekin, ne, K1, K2, K3, X))) {
      std::cout << "   new low is " << newVal << std::endl;
      low = newVal;
    }
    else {
      std::cout << "   new high is " << newVal << std::endl;
      high = newVal;
    }
    iter++;
  }
  bm = newVal;
  bn = newVal;

  //// use beta to find chemical potential
  double mue = 0.0;
  double nue = 4*ne*pow(M_PI*bm/(2*p->me),1.5);	// constant to simplify
  mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bm;
  std::cout << "Chemical potential " << mue << std::endl;

  // TODO account for temperature dropping in time
  std::cout << "inverse temp is " << bn << std::endl;
  std::cout << std::endl;
  for (int ii = 0; ii < p->Nk; ii++) {
    // TODO factor in Boltzmann constant?
    fdd[ii] = 1.0/(1.0 + exp((ekin - mue)*bn));
    std::cout << 1.0/(1.0 + exp((ekin - mue)*bn)) << " ";
    std::cout << "FDD[" << ii << "]: " << std::scientific << fdd[ii] << std::endl;
  }

  return;
}

/* implements equation B13 from Binder et. al, PRB 1991.
 * bm is the beta (1/kT) value.
 * ekin is the kinetic energy
 * ne is the carrier density
 * K1-3 are constants
 * X is a constant
 */
double b13(double bm, double ekin, double ne, double K1, double K2, double K3, double X) {
  return -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
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
  std::cout << "POPULATION " << NV_Ith_S(y, 0) << std::endl;
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

