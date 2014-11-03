#include "dynamix.hpp"

// DEBUG compiler flag: turn on to generate basic debug outputs.
// #define DEBUG

// DEBUG2 flag: turn on for more numerical output
// #define DEBUG2

// #define DEBUG_UPDATE

/* Updates \rho(t) at each time step. */
void updateDM(N_Vector dm, int timeStep, Params * p) {
  timeStep = timeStep*p->numOutputSteps/p->numsteps;
#ifdef DEBUG_UPDATE
  std::cout << "Updating DM at step " << timeStep << "...";
#endif
  int N = 2*p->NEQ2;
  memcpy(&(p->dmt[N*timeStep]), N_VGetArrayPointer(dm), N*sizeof(realtype));
#ifdef DEBUG_UPDATE
  std::cout << "done.\n";
#endif

  return;
}

/* Updates \psi(t) at each time step. */
void updateWfn(N_Vector wfn, int timeStep, Params * p) {
  timeStep = timeStep*p->numOutputSteps/p->numsteps;
#ifdef DEBUG_UPDATE
  std::cout << "Updating wavefunction at time step " << timeStep << "..." << std::endl;
  std::cout << "Wavefunction is " << std::endl;
  N_VPrint_Serial(wfn);
#endif
  int N = 2*p->NEQ;
  memcpy(&(p->wfnt[N*timeStep]), N_VGetArrayPointer(wfn), N*sizeof(realtype));
#ifdef DEBUG_UPDATE
  std::cout << "done updating wavefunction." << std::endl;
#endif
  return;
}

/* Get band index based on flag */
int bandStartIdx(int bandFlag, Params * p) {
  if (bandFlag == CONDUCTION) {
    return p->Ik;
  }
  else if (bandFlag == VALENCE) {
    return p->Il;
  }
  else if (bandFlag == BRIDGE) {
    return p->Ib;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Ic;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}

/* Get band index based on flag */
int bandEndIdx(int bandFlag, Params * p) {
  if (bandFlag == CONDUCTION) {
    return p->Ik + p->Nk;
  }
  else if (bandFlag == VALENCE) {
    return p->Il + p->Nl;
  }
  else if (bandFlag == BRIDGE) {
    return p->Ib + p->Nb;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Ic + p->Nc;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}

/* get number of states in band */
int bandNumStates(int bandFlag, Params * p) {
  if (bandFlag == CONDUCTION) {
    return p->Nk;
  }
  else if (bandFlag == VALENCE) {
    return p->Nl;
  }
  else if (bandFlag == BRIDGE) {
    return p->Nb;
  }
  else if (bandFlag == QD_CONDUCTION) {
    return p->Nc;
  }
  else {
    std::cout << "WARNING [" << __FUNCTION__ << "]: unspecified band with flag " << bandFlag << std::endl;
    return 0;
  }
}

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, Params * p) {
  int s;  // sign +/-
  double m; // mass of electron/hole

  // determine conduction vs. valence band, electron/hole masses
  if (flag == CONDUCTION) {
    s = 1;
    m = p->me;
  }
  else { // (flag == VALENCE)
    s = -1;
    m = p->mh;
  }

  // assign energies
  for (int ii = 0; ii < n; ii++) {
    energies[ii] = bandEdge + s*ii*ii/(2*m*pow(p->X2,2));
  }

  return;
}

/* Updates the Hamiltonian with the time-dependent torsional coupling
 * and laser field.
 */
void updateHamiltonian(Params * p, realtype t) {
  // TODO unpack NEQ from p

  // get pointer to H
  realtype * H = &(p->H)[0];

  //// first handle torsion
  double torsionValue = 0.0;
  if (p->torsion) {
    torsionValue = p->getTorsionCoupling(t);
#ifdef DEBUG_TORSION
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

#ifdef __USE_MKL__
  // make sparse version of Hamiltonian
  int job [6] = {0, 0, 0, 2, p->NEQ2, 1};
  int info = 0;

  mkl_ddnscsr(&job[0], &(p->NEQ), &(p->NEQ), &(p->H)[0], &(p->NEQ), &(p->H_sp)[0],
      &(p->H_cols)[0], &(p->H_rowind)[0], &info);
#endif

  return;
}

void initialize(Params * p, const bool readFiles) {
  initHamiltonian(p, readFiles);
  initWavefunction(p);
}

// This function builds up the Hamiltonian, as well as the constituent site
// energy array and coupling array.
void initHamiltonian(Params * p, const bool readFiles) {
  if (readFiles) {
    // first read in energies and couplings from files
    p->Nc = numberOfValuesInFile(p->cEnergiesInput.c_str());
    p->Nb = numberOfValuesInFile(p->bEnergiesInput.c_str());

    p->k_energies.resize(p->Nk, 0.0);
    p->c_energies.resize(p->Nc, 0.0);
    p->b_energies.resize(p->Nb, 0.0);
    p->l_energies.resize(p->Nl, 0.0);

    // c energies are defined in file
    readVectorFromFile(p->c_energies, p->cEnergiesInput.c_str(), p->Nc);

    // bridge-dependent parameters
    if (p->bridge_on) {
      if (p->Nb < 1) {
        std::cerr << "\nERROR: bridge_on but no bridge states.  The file b_energies.in is probably empty.\n";
        exit(-1);
      }

      readVectorFromFile(p->b_energies, p->bEnergiesInput.c_str(), p->Nb);

      p->Vbridge.resize(p->Nb+1);
      readVectorFromFile(p->Vbridge, p->VBridgeInput.c_str(), p->Nb + 1);

  #ifdef DEBUG
      std::cout << "COUPLINGS:";
      for (int ii = 0; ii < p->Nb+1; ii++) {
        std::cout << " " << p->Vbridge[ii];
      }
      std::cout << std::endl;
  #endif
    }
    else { // no bridge
      p->Nb = 0;
      p->Vnobridge.resize(1);
      readVectorFromFile(p->Vnobridge, p->VNoBridgeInput.c_str(), 1);
    }
  }

#ifdef DEBUG
  std::cout << "\nDone reading things from inputs.\n";
#endif

  // assign bulk conduction and valence band energies
  buildContinuum(&(p->k_energies[0]), p->Nk, p->kBandEdge, p->kBandTop);
  buildContinuum(&(p->l_energies[0]), p->Nl, p->kBandEdge - p->valenceBand - p->bulk_gap, p->kBandEdge - p->bulk_gap);

  // calculate band width
  p->kBandWidth = fabs(p->k_energies.back() - p->k_energies.front());

  // set total number of equations
  p->NEQ = p->Nk+p->Nc+p->Nb+p->Nl;                          // total number of equations set
  p->NEQ2 = p->NEQ*p->NEQ;                         // number of elements in DM

  // set index start positions for each type of state
  p->Ik = 0;
  p->Ic = p->Nk;
  p->Ib = p->Ic+p->Nc;
  p->Il = p->Ib+p->Nb;


  //// ASSEMBLE ARRAY OF ENERGIES


  p->energies.resize(p->NEQ);
  for (int ii = 0; ii < p->Nk; ii++) {
    p->energies[p->Ik + ii] = p->k_energies[ii];
  }
  for (int ii = 0; ii < p->Nc; ii++) {
    p->energies[p->Ic + ii] = p->c_energies[ii];
  }
  for (int ii = 0; ii < p->Nb; ii++) {
    p->energies[p->Ib + ii] = p->b_energies[ii];
  }
  for (int ii = 0; ii < p->Nl; ii++) {
    p->energies[p->Il + ii] = p->l_energies[ii];
  }

#ifdef DEBUG
  for (int ii = 0; ii < p->NEQ; ii++) {
    std::cout << "energies[" << ii << "] is " << p->energies[ii] << "\n";
  }
#endif

  // set up Hamiltonian ////////////////////////////////////////////////////////

#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif

  p->buildCoupling();
  p->buildHamiltonian();
}

// This function builds up the initial wavefunction coefficients based on inputs
void initWavefunction(Params * p) {
  std::vector<realtype> k_coeffs (p->Nk, 0.0);
  std::vector<realtype> c_coeffs (p->Nc, 0.0);
  std::vector<realtype> b_coeffs (p->Nb, 0.0);
  std::vector<realtype> l_coeffs (p->Nl, 0.0);

  double summ;


  //// BUILD INITIAL WAVEFUNCTION

  // set coefficients in each band of states

  // bulk valence band /////////////////////////////////////////////////////////
  if (p->VBPopFlag == POP_EMPTY) {
#ifdef DEBUG
    std::cout << "Initializing empty valence band" << std::endl;
#endif
    l_coeffs.assign(l_coeffs.size(), 0.0);
  }
  else if (p->VBPopFlag == POP_FULL) {
#ifdef DEBUG
    std::cout << "Initializing full valence band" << std::endl;
#endif
    l_coeffs.assign(l_coeffs.size(), 1.0);
  }
  else {
    std::cerr << "ERROR: unrecognized VBPopFlag " << p->VBPopFlag << std::endl;
  }

  // bulk conduction band //////////////////////////////////////////////////////
  if (p->CBPopFlag == POP_EMPTY) {
#ifdef DEBUG
    std::cout << "Initializing empty conduction band" << std::endl;
#endif
    k_coeffs.assign(k_coeffs.size(), 0.0);
  }
  else if (p->CBPopFlag == POP_FULL) {
#ifdef DEBUG
    std::cout << "Initializing full conduction band" << std::endl;
#endif
    k_coeffs.assign(k_coeffs.size(), 1.0);
  }
  else if (p->CBPopFlag == POP_CONSTANT) {
#ifdef DEBUG
    std::cout << "Initializing constant distribution in conduction band" << std::endl;
#endif
    k_coeffs.assign(k_coeffs.size(), 0.0);
    initializeArray(&(k_coeffs[p->Nk_first-1]), p->Nk_final - p->Nk_first + 1, 1.0);
  }
  else if (p->CBPopFlag == POP_GAUSSIAN) {
#ifdef DEBUG
    std::cout << "Initializing Gaussian in conduction band" << std::endl;
#endif
    buildKPopsGaussian(&(k_coeffs[0]), &(p->energies[p->Ik]), p->kBandEdge, p->bulkGaussSigma, p->bulkGaussMu, p->Nk);
  }
  else {
    std::cerr << "ERROR: unrecognized CBPopFlag " << p->CBPopFlag << std::endl;
  }

  //// QD //////////////////////////////////////////////////////////////////////
  if (p->QDPopFlag == POP_EMPTY) {
    c_coeffs.assign(c_coeffs.size(), 0.0);
  }
  else if (p->QDPopFlag == POP_FULL) {
    c_coeffs.assign(c_coeffs.size(), 1.0);
  }
  else if (p->QDPopFlag == POP_CONSTANT) {
#ifdef DEBUG
    std::cout << "Initializing constant distribution in QD band" << std::endl;
#endif
    c_coeffs.assign(c_coeffs.size(), 0.0);
    initializeArray(&(c_coeffs[p->Nc_first-1]), p->Nc_final - p->Nc_first + 1, 1.0);
  }
  else {
    std::cerr << "ERROR: unrecognized QDPopFlag " << p->QDPopFlag << std::endl;
  }

  // create empty wavefunction
  // resize method does not initialize all values, hence second call
  p->startWfn.resize(2*p->NEQ, 0.0);
  initializeArray(&(p->startWfn[0]), p->startWfn.size(), 0.0);
  if (!p->wavefunction) {
    p->startDM.resize(2*p->NEQ2, 0.0);
    initializeArray(&(p->startDM[0]), p->startDM.size(), 0.0);
  }

  // assign real parts of wavefunction coefficients (imaginary are zero)
  for (int ii = 0; ii < p->Nk; ii++) {
    p->startWfn[p->Ik + ii] = k_coeffs[ii];
  }
  for (int ii = 0; ii < p->Nc; ii++) {
    p->startWfn[p->Ic + ii] = c_coeffs[ii];
  }
  for (int ii = 0; ii < p->Nb; ii++) {
    p->startWfn[p->Ib + ii] = b_coeffs[ii];
  }
  for (int ii = 0; ii < p->Nl; ii++) {
    p->startWfn[p->Il + ii] = l_coeffs[ii];
  }

  if (isOutput(p->outs, "psi_start.out")) {
    outputWavefunction(&(p->startWfn[0]), p->NEQ);
  }

  // Give all coefficients a random phase
  if (p->random_phase) {
    float phi;
    // set the seed
    if (p->random_seed == -1) { srand(time(NULL)); }
    else { srand(p->random_seed); }
    for (int ii = 0; ii < p->NEQ; ii++) {
      phi = 2*3.1415926535*(float)rand()/(float)RAND_MAX;
      p->startWfn[ii] = p->startWfn[ii]*cos(phi);
      p->startWfn[ii + p->NEQ] = p->startWfn[ii + p->NEQ]*sin(phi);
    }
  }

#ifdef DEBUG
  // print out details of wavefunction coefficients
  std::cout << std::endl;
  for (int ii = 0; ii < p->Nk; ii++) {
    std::cout << "starting wavefunction: Re[k(" << ii << ")] = " << p->startWfn[p->Ik + ii] << std::endl;
  }
  for (int ii = 0; ii < p->Nc; ii++) {
    std::cout << "starting wavefunction: Re[c(" << ii << ")] = " << p->startWfn[p->Ic + ii] << std::endl;
  }
  for (int ii = 0; ii < p->Nb; ii++) {
    std::cout << "starting wavefunction: Re[b(" << ii << ")] = " << p->startWfn[p->Ib + ii] << std::endl;
  }
  for (int ii = 0; ii < p->Nl; ii++) {
    std::cout << "starting wavefunction: Re[l(" << ii << ")] = " << p->startWfn[p->Il + ii] << std::endl;
  }
  for (int ii = 0; ii < p->Nk; ii++) {
    std::cout << "starting wavefunction: Im[k(" << ii << ")] = " << p->startWfn[p->Ik + ii + p->NEQ] << std::endl;
  }
  for (int ii = 0; ii < p->Nc; ii++) {
    std::cout << "starting wavefunction: Im[c(" << ii << ")] = " << p->startWfn[p->Ic + ii + p->NEQ] << std::endl;
  }
  for (int ii = 0; ii < p->Nb; ii++) {
    std::cout << "starting wavefunction: Im[b(" << ii << ")] = " << p->startWfn[p->Ib + ii + p->NEQ] << std::endl;
  }
  for (int ii = 0; ii < p->Nl; ii++) {
    std::cout << "starting wavefunction: Im[l(" << ii << ")] = " << p->startWfn[p->Il + ii + p->NEQ] << std::endl;
  }
  std::cout << std::endl;
  summ = 0;
  for (int ii = 0; ii < 2*p->NEQ; ii++) {
    summ += pow(p->startWfn[ii],2);
  }
  std::cout << "\nTotal population is " << summ << "\n\n";
#endif

  //// CREATE DENSITY MATRIX

  if (! p->wavefunction) {
    initializeArray(&(p->startDM[0]), p->startDM.size(), 0.0);
#pragma omp parallel for
    for (int ii = 0; ii < p->NEQ; ii++) {
      // diagonal part
      p->startDM[p->NEQ*ii + ii] = pow(p->startWfn[ii],2) + pow(p->startWfn[ii + p->NEQ],2);
      if (p->coherent) {
#ifdef DEBUG
        if (ii == 0) {
          std::cout << "\nDM starting state is coherent.\n\n";
        }
#endif
        // off-diagonal part
        for (int jj = 0; jj < ii; jj++) {
          // real part of \rho_{ii,jj}
          p->startDM[p->NEQ*ii + jj] = p->startWfn[ii]*p->startWfn[jj] + p->startWfn[ii+p->NEQ]*p->startWfn[jj+p->NEQ];
          // imaginary part of \rho_{ii,jj}
          p->startDM[p->NEQ*ii + jj + p->NEQ2] = p->startWfn[ii]*p->startWfn[jj+p->NEQ] - p->startWfn[jj]*p->startWfn[ii+p->NEQ];
          // real part of \rho_{jj,ii}
          p->startDM[p->NEQ*jj + ii] = p->startDM[p->NEQ*ii + jj];
          // imaginary part of \rho_{jj,ii}
          p->startDM[p->NEQ*jj + ii + p->NEQ2] = -1*p->startDM[p->NEQ*ii + jj + p->NEQ*p->NEQ];
        }
      }
      else {
#ifdef DEBUG
        if (ii == 0) {
          std::cout << "\nDM starting state is incoherent.\n\n";
        }
#endif
      }
    }

#ifdef DEBUG2
    // print out density matrix
    std::cout << "\nDensity matrix without normalization:\n\n";
    for (int ii = 0; ii < p->NEQ; ii++) {
      for (int jj = 0; jj < p->NEQ; jj++) {
        fprintf(stdout, "(%+.1e,%+.1e) ", p->startDM[p->NEQ*ii + jj], p->startDM[p->NEQ*ii + jj + p->NEQ2]);
      }
      fprintf(stdout, "\n");
    }
#endif

    // Normalize the DM so that populations add up to 1.
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      // assume here that diagonal elements are all real
      summ += p->startDM[p->NEQ*ii + ii];
    }
    if ( summ == 0.0 ) {
      std::cerr << "\nFATAL ERROR [populations]: total population is 0!\n";
      exit(-1);
    }
    if (summ != 1.0) {
      // the variable 'summ' is now a multiplicative normalization factor
      summ = 1.0/summ;
      for (int ii = 0; ii < 2*p->NEQ2; ii++) {
        p->startDM[ii] *= summ;
      }
    }
#ifdef DEBUG
    std::cout << "\nThe normalization factor for the density matrix is " << summ << "\n\n";
#endif

    // Error checking for total population; recount population first
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += p->startDM[p->NEQ*ii + ii];
    }
    if ( fabs(summ-1.0) > 1e-12) {
      std::cerr << "\nWARNING [populations]: After normalization, total population is not 1, it is " << summ << "!\n";
    }
#ifdef DEBUG
    std::cout << "\nAfter normalization, the sum of the populations in the density matrix is " << summ << "\n\n";
#endif
  }
  // wavefunction
  else {
    // normalize
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->startWfn[ii],2) + pow(p->startWfn[ii+p->NEQ],2);
    }
#ifdef DEBUG
    std::cout << "Before normalization, the total population is " << summ << std::endl;
#endif
    summ = 1.0/sqrt(summ);
    for (int ii = 0; ii < 2*p->NEQ; ii++) {
      p->startWfn[ii] *= summ;
    }

    // check total population
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += pow(p->startWfn[ii],2) + pow(p->startWfn[ii+p->NEQ],2);
    }
#ifdef DEBUG
    std::cout << "After normalization, the total population is " << summ << std::endl;
#endif
    if (fabs(summ - 1.0) > 1e-12) {
      std::cerr << "WARNING: wavefunction not normalized!  Total density is " << summ << std::endl;
    }
  }
}
