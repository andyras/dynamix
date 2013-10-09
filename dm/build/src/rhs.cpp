#include "rhs.hpp"

//#define DEBUG_RHS
//#define DEBUG_RTA
//
// DEBUGf flag: general output at each CVode step
//#define DEBUGf
//
// DEBUGf_DM flag: DEBUGf for density matrix EOM
//#define DEBUGf_DM


/* Updates the Hamiltonian with the time-dependent torsional coupling
 * and laser field.
 */
void updateHamiltonian(PARAMETERS * p, realtype t) {
  // get pointer to H
  realtype * H = &(p->H)[0];

  //// first handle torsion
  if (p->torsion) {
    double torsionValue = p->torsionV->value(t);
#ifdef DEBUG_RHS
    std::cout << "Value of torsion-mediated coupling is " << torsionValue << std::endl;
#endif

    // bridge is off, coupling is between k and c states
    if (!(p->bridge_on)) {
#ifdef DEBUG_RHS
      std::cout << "torsion between k and c states" << std::endl;
#endif
      for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
	for (int jj = p->Ic; jj < (p->Ic + p->Nc); jj++) {
	  H[ii*p->NEQ + jj] = torsionValue;
	  H[jj*p->NEQ + ii] = torsionValue;
	}
      }
    }
    // torsion is at first bridge coupling
    else if (p->torsionSite == 0) {
#ifdef DEBUG_RHS
      std::cout << "torsion between k states and bridge" << std::endl;
#endif
      for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
	H[ii*p->NEQ + p->Ib] = torsionValue;
	H[p->Ib*p->NEQ + ii] = torsionValue;
      }
    }
    // torsion is at last bridge coupling
    else if (p->torsionSite == p->Nb) {
#ifdef DEBUG_RHS
      std::cout << "torsion between bridge and c states" << std::endl;
#endif
      for (int ii = p->Ic; ii < (p->Ic + p->Nc); ii++) {
	H[ii*p->NEQ + p->Ib + p->Nb - 1] = torsionValue;
	H[(p->Ib + p->Nb - 1)*p->NEQ + ii] = torsionValue;
      }
    }
    // torsion is between bridge sites
    else {
#ifdef DEBUG_RHS
      std::cout << "torsion between bridge sites " << p->torsionSite - 1
	<< " and " << p->torsionSite << "." << std::endl;
#endif
      H[(p->Ib + p->torsionSite - 1)*p->NEQ + p->Ib + p->torsionSite] = torsionValue;
      H[(p->Ib + p->torsionSite)*p->NEQ + p->Ib + p->torsionSite - 1] = torsionValue;
    }
  }

  //// now handle pump pulse
  double laserCoupling = 0.0;
  if (p->laser_on) {
    laserCoupling = gaussPulse(t, p->pumpFWHM, p->pumpAmpl, p->pumpPeak, p->pumpFreq, p->pumpPhase);
#ifdef DEBUG_RHS
    std::cout << "Value of laser coupling between valence and conduction bands is " << laserCoupling << std::endl;
#endif
    // coupling is between valence and conduction bands
    for (int ii = p->Il; ii < (p->Il + p->Nl); ii++) {
      for (int jj = p->Ik; jj < (p->Ik + p->Nk); jj++) {
	H[(ii)*p->NEQ + jj] = laserCoupling;
	H[(jj)*p->NEQ + ii] = laserCoupling;
      }
    }
  }

  return;
}

/* Right-hand-side equation for density matrix */
int RHS_DM(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;

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

  std::cout << "Closing output file for density matrix coefficients in time.\n";
  fclose(dmf);
#endif

  return 0;
}

/* gives the equilibrated FDD for the system */
void buildFDD(struct PARAMETERS * p, N_Vector y, double * fdd) {
  //// "fine structure constant" -- conversion from index to wave vector
#ifdef DEBUG_RTA
  std::cout << "p->X2   " << p->X2 << std::endl;
#endif

  //// calculate n_e and e_kin
  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(p->X2,3));
#ifdef DEBUG_RTA
  std::cout << "factor   " << factor << std::endl;
#endif
  // assign vector of energies
  // std::vector<double> E (p->Nk,0.0);
  realtype * E = new realtype [p->Nk];
  for (int ii = 0; ii < p->Nk; ii++) {
    E[ii] = pow(ii,2)/(2*p->me*pow(p->X2,2));
  }

#ifdef DEBUG_RTA
  std::cout << std::setprecision(28);
#endif
  // Simpson's Rule method
  // skip the first point because the value will be zero
  double SF = 4.0;	// Simpson's factor
  int sign = -1;	// sign
  for (int ii = 1; ii < (p->Nk-1); ii++) {
    ne += SF*factor*ii*ii*NV_Ith_S(y, ii*p->NEQ + ii);
    ekin += SF*factor*pow(ii,2)*NV_Ith_S(y, ii*p->NEQ + ii)*E[ii];
#ifdef DEBUG_RTA
    std::cout << "Ne " << ii*ii << "*" << SF << "/3.0*" << NV_Ith_S(y, ii*p->NEQ + ii) << "/" << pow(5.29e-11,3)/factor << std::endl;
    std::cout << "ekin " << pow(ii,4) << "*" << SF << "/3.0*" << NV_Ith_S(y, ii*p->NEQ + ii) << "*" << 4.3597482e-18/(2*p->me*pow(p->X2,2)) << "/" << pow(5.29e-11,3)/factor << std::endl;
    std::cout << "ekin " << ekin/pow(5.29e-11,3)*4.3597482e-18/3.0 
      << " += " << SF*factor*pow(ii,4)*NV_Ith_S(y, ii*p->NEQ + ii)/(2*p->me*pow(p->X2,2))/pow(5.29e-11,3)*4.3597482e-18/3.0 << std::endl;
#endif
    SF += sign*2.0;
    sign *= -1;
  }
  // add last point
  ne += factor*pow(p->Nk-1,2)*NV_Ith_S(y, (p->Nk - 1)*p->NEQ + p->Nk - 1);
  ekin += factor*pow(p->Nk-1,2)*NV_Ith_S(y, (p->Nk - 1)*p->NEQ + p->Nk - 1)*E[p->Nk - 1];
  // divide by three
  ne /= 3.0;
  ekin /= 3.0;

#ifdef DEBUG_RTA
  std::cout << "ne        " << ne << std::endl;
  std::cout << "ne (SI)   " << ne/pow(5.29e-11,3) << std::endl;
  std::cout << "ekin      " << ekin << std::endl;
  std::cout << "ekin (SI) " << ekin/pow(5.29e-11,3)*4.3597482e-18 << std::endl;
#endif

  //// find the inverse temperature (beta)
  int iter = 0;
  const int maxiter = 60;
  double tol = 1e-12;
  double K1 = 4.8966851;		// constants
  double K2 = 0.04496457;
  double K3 = 0.133376;
  double X = 4*ne*pow(M_PI/(2*p->me),1.5)*6.9608/6.95369; // FIXME conversion at end to match Sai's values...
#ifdef DEBUG_RTA
  std::cout << "XX " << X/pow(2.293710449e+17,1.5) << std::endl;
#endif
  double bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
  double bm = 0.0;
  double vol = pow(1.0/5.29e-11,3);		// volume element, multiply to go from a0^-3 to m^-3

  // loop applies Newton-Raphson method to get zero of function
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
#ifdef DEBUG_RTA
  std::cout << "Newton-Raphson to find inverse temperature" << std::endl;
#endif
  while ((fabs(bn - bm)/bm > tol) && (iter < maxiter)) {
    bm = bn;
    f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
    fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*X*pow(bm,0.5));
#ifdef DEBUG_RTA
    std::cout << "Iteration     " << std::setw(15) << iter << std::endl;
    std::cout << "bm            " << std::setw(15) << bm/4.3597482e-18 << std::endl;
    std::cout << "f(bm) term 1: " << std::setw(15) << vol*-bm*ekin << std::endl;
    std::cout << "f(bm) term 2: " << std::setw(15) << vol*1.5*ne*(1 + K1) << std::endl;
    std::cout << "f(bm) term 3: " << std::setw(15) << vol*1.5*ne*(-1*K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5))) << std::endl;
    std::cout << "f(bm) term 4: " << std::setw(15) << vol*1.5*ne*(0.5*K3*X*pow(bm,1.5)) << std::endl;
    std::cout << "f(bm) (SI)    " << std::setw(15) << f*pow(1.0/5.29e-11,3) << std::endl;
    std::cout << "f(bm) (a.u)   " << std::setw(15) << f << std::endl;
    std::cout << "f'(bm) (SI)   " << std::setw(15) << fp*pow(1.0/5.29e-11,3)*4.3597482e-18 << std::endl;
    std::cout << "f'(bm) (a.u)  " << std::setw(15) << fp << std::endl;
#endif
    bn = bm - f/fp;
    iter++;
  }
#ifdef DEBUG_RTA
  std::cout << std::endl;
#endif

  //// use beta to find chemical potential
  double mue = 0.0;
  double nue = 4*ne*pow(M_PI*bn/(2*p->me),1.5);	// constant to simplify
  mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bn;
#ifdef DEBUG_RTA
  std::cout << "Chemical potential " << mue*4.3597482e-18 << std::endl;
#endif

  // TODO account for temperature dropping in time
#ifdef DEBUG_RTA
  std::cout << "inverse temp is " << bn << std::endl;
  std::cout << std::endl;
#endif

  //// assign Fermi-Dirac function
  for (int ii = 0; ii < p->Nk; ii++) {
    // TODO factor in Boltzmann constant?
    fdd[ii] = 1.0/(1.0 + exp((E[ii] - mue)*bn));
#ifdef DEBUG_RTA
    std::cout << "FDD[" << ii << "]: " << std::scientific << fdd[ii] << std::endl;
#endif
  }

  // free array
  delete [] E;

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

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif

#ifdef DEBUG_RHS
  std::cout << "Time " << t << std::endl;
#endif

  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  //std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  realtype * H = &(p->H)[0];
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g1 = p->gamma1;
  realtype g2 = p->gamma2;

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
    NV_Ith_S(ydot, ii) = 0.0;
  }

  //// diagonal; no need to calculate the imaginary part
  //   get equilibrium FDD populations
  //std::vector<double> fdd(p->Nk);
  double * fdd = new double [p->Nk];
#ifdef DEBUG_RTA
  std::cout << "POPULATION " << NV_Ith_S(y, 0) << std::endl;
#endif
  buildFDD(p, y, fdd);

  //// normalize FDD to amount of population in conduction band
  double fddSum = 0.0;
  double CBSum = 0.0;
  for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
    // sum population in FDD
    fddSum += fdd[ii - p->Ik];
    // sum population in CB
    CBSum += NV_Ith_S(y, ii*N + ii);
  }
  double fddNorm = CBSum/fddSum;
#ifdef DEBUG_RTA
  std::cout << "FDD normalization constant is " << fddNorm << std::endl;
#endif
  for (int ii = 0; ii < p->Nk; ii++) {
    fdd[ii] *= fddNorm;
  }

#pragma omp parallel for
  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      NV_Ith_S(ydot, ii*N + ii) += 2*H[ii*N + jj]*NV_Ith_S(y, jj*N + ii + N2);
    }
  }
  // force conduction band toward Fermi-Dirac distribution
  for (int ii = p->Ik; ii < (p->Ik + p->Nk); ii++) {
    NV_Ith_S(ydot, ii*N + ii) -= g1*(NV_Ith_S(y, ii*N + ii) - fdd[ii]);
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

  // free fdd
  delete [] fdd;

  return 0;
}

/* Right-hand-side equation for density matrix
 * using dephasing */
int RHS_DM_dephasing(realtype t, N_Vector y, N_Vector ydot, void * data) {

#ifdef DEBUGf_DM
  // file for density matrix coeff derivatives in time
  FILE * dmf;
  std::cout << "Creating output file for density matrix coefficient derivatives in time.\n";
  dmf = fopen("dmf.out", "w");
#endif


  // data is a pointer to the params struct
  PARAMETERS * p;
  p = (PARAMETERS *) data;

  // extract parameters from p
  std::vector<realtype> H = p->H; // copying vector is OK performance-wise
  int N = p->NEQ;
  int N2 = p->NEQ2;
  realtype g2 = p->gamma2;

  // update Hamiltonian if it is time-dependent
  if (p->torsion || p->laser_on) {
    updateHamiltonian(p, t);
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

