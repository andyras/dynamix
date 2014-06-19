#include "dynamix.hpp"

// DEBUG compiler flag: turn on to generate basic debug outputs.
#define DEBUG

// DEBUG2 flag: turn on for more numerical output
// #define DEBUG2

// DEBUGf: debug inner CVode loop
// #define DEBUGf

#define DEBUG_BUILDCOUPLING
#define DEBUG_UPDATEDM
#define DEBUG_UPDATEWFN
//#define DEBUG_BUILDHAMILTONIAN

/* builds a Hamiltonian from site energies and couplings. */
void buildHamiltonian(realtype * H, std::vector<realtype> & energy, std::vector< std::vector<realtype> > * V, struct PARAMETERS * p) {
  // indices
  int idx1, idx2;
  int N = p->NEQ;

#ifdef DEBUG_BUILDHAMILTONIAN
  fprintf(stderr, "Assigning diagonal elements of Hamiltonian.\n");
#endif
  for (int ii = 0; ii < N; ii++) {
    // diagonal
    H[ii*N + ii] = energy[ii];
#ifdef DEBUG_BUILDHAMILTONIAN
    std::cout << "diagonal element " << ii << " of H is " << energy[ii] << "\n";
#endif
  }

  if (p->bridge_on) {
    // assign bulk-bridge coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-bridge coupling elements in Hamiltonian.\n");
#endif
    idx2 = p->Ib;
    for (int ii = 0; ii < p->Nk; ii++) {
      idx1 = p->Ik + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V->at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V->at(idx1)[idx2];
      H[idx2*N + idx1] = V->at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bulk-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-bridge couplings
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-bridge coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < (p->Nb-1); ii++) {
      idx1 = p->Ib + ii;
      idx2 = p->Ib+ ii + 1;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V->at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V->at(idx1)[idx2];
      H[idx2*N + idx1] = V->at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bridge-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-QD coupling elements in Hamiltonian.\n");
#endif
    idx2 = p->Ib + p->Nb - 1;
    for (int ii = 0; ii < p->Nc; ii++) {
      idx1 = p->Ic + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V->at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V->at(idx1)[idx2];
      H[idx2*N + idx1] = V->at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bridge-QD coupling elements in Hamiltonian.\n");
  }
  // no bridge
  else {
    // assign bulk-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-QD coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < p->Nk; ii++) {
      idx1 = p->Ik + ii;
      for (int jj = 0; jj < p->Nc; jj++) {
  idx2 = p->Ic + jj;
#ifdef DEBUG_BUILDHAMILTONIAN
  fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
  fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
  fprintf(stderr, "%e\n", V->at(idx1)[idx2]);
#endif
  H[idx1*N + idx2] = V->at(idx1)[idx2];
  H[idx2*N + idx1] = V->at(idx2)[idx1];
      }
    }
  }
}

/* Updates \rho(t) at each time step. */
void updateDM(N_Vector dm, realtype * dmt, int timeStep, struct PARAMETERS * p) {
#ifdef DEBUG_UPDATEDM
  std::cout << "Updating DM at step " << timeStep << "...";
#endif
  int N = 2*p->NEQ2;
  memcpy(&dmt[N*timeStep], N_VGetArrayPointer(dm), N*sizeof(realtype));
  /*
     realtype * nv = N_VGetArrayPointer(dm);
     nv_tp = &dmt[N*timeStep];
     for (int ii = 0; ii < N; ii++) {
     nv_tp[ii] = nv[ii];
     }
     */
#ifdef DEBUG_UPDATEDM
  std::cout << "done.\n";
#endif

  return;
}

/* Updates \psi(t) at each time step. */
void updateWfn(N_Vector wfn, realtype * wfnt, int timeStep, struct PARAMETERS * p) {
#ifdef DEBUG_UPDATEWFN
  std::cout << "Updating wavefunction at time step " << timeStep << "..." << std::endl;
  std::cout << "Wavefunction is " << std::endl;
  N_VPrint_Serial(wfn);
#endif
  int N = 2*p->NEQ;
  memcpy(&wfnt[N*timeStep], N_VGetArrayPointer(wfn), N*sizeof(realtype));
#ifdef DEBUG_UPDATEWFN
  std::cout << "done updating wavefunction." << std::endl;
#endif
  return;
}

/* assign coupling constants to global array V */
void buildCoupling (std::vector< std::vector<realtype> > * vArray, struct PARAMETERS * p, std::map<const std::string, bool> &outs) {

  double Vkc; // coupling between bulk and QD
  double Vkb1;  // coupling between bulk and first bridge
  double VbNc;  // coupling between last bridge and QD

  // initialize the coupling array
  vArray->resize(p->NEQ);
  for (int ii = 0; ii < vArray->size(); ii++) {
    vArray->at(ii).assign(p->NEQ, 0.0);
  }

std::cout << "\n\n\nWHOOOOOT\n\n\n";

#ifdef DEBUG_BUILDCOUPLING
  if (p->bridge_on) {
    for (int ii = 0; ii < p->Nb + 1; ii++) {
      std::cout << "p->Vbridge[" << ii << "] is ";
      std::cout << p->Vbridge[ii] << "\n";
    }
  }
#endif

  // bridge
  if (p->bridge_on) {
    // coupling between k and b1
    if ((p->scale_bubr) && (p->Nk > 1)) {
      Vkb1 = sqrt(p->Vbridge[0]*(p->kBandTop-p->kBandEdge)/(p->Nk-1));
    }
    else {
      Vkb1 = p->Vbridge[0];
    }
    if (p->parabolicCoupling) {
      for (int ii = 0; ii < p->Nk; ii++) {
        vArray->at(p->Ik+ii)[p->Ib] = parabolicV(Vkb1, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
        vArray->at(p->Ib)[p->Ik+ii] = parabolicV(Vkb1, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
      }
    }
    else {
      for (int ii = 0; ii < p->Nk; ii++) {
        vArray->at(p->Ik+ii)[p->Ib] = Vkb1;
        vArray->at(p->Ib)[p->Ik+ii] = Vkb1;
      }
    }

    // coupling between bN and c
    if ((p->scale_brqd) && (p->Nc > 1)) {
      VbNc = p->Vbridge[p->Nb]/sqrt(p->Nc-1);
    }
    else {
      VbNc = p->Vbridge[p->Nb];
    }
    for (int ii = 0; ii < p->Nc; ii++) {
      vArray->at(p->Ic+ii)[p->Ib+p->Nb-1] = VbNc;
      vArray->at(p->Ib+p->Nb-1)[p->Ic+ii] = VbNc;
    }

    // coupling between bridge states
    for (int ii = 0; ii < p->Nb - 1; ii++) {
      vArray->at(p->Ib+ii)[p->Ib+ii+1] = p->Vbridge[ii+1];
      vArray->at(p->Ib+ii+1)[p->Ib+ii] = p->Vbridge[ii+1];
    }
  }
  // no bridge
  else {
    // scaling
    if ((p->scale_buqd) && (p->Nk > 1)) {
      Vkc = sqrt(p->Vnobridge[0]*(p->kBandTop-p->kBandEdge)/(p->Nk-1));
    }
    else {
      Vkc = p->Vnobridge[0];
    }

    // parabolic coupling of bulk band to QD
    if (p->parabolicCoupling) {
      for (int ii = 0; ii < p->Nk; ii++) {
  for (int jj = 0; jj < p->Nc; jj++) {
    vArray->at(p->Ik+ii)[p->Ic+jj] = parabolicV(Vkc, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
    vArray->at(p->Ic+jj)[p->Ik+ii] = parabolicV(Vkc, p->energies[p->Ik+ii], p->kBandEdge, p->kBandTop);
  }
      }
    }
    else {
      for (int ii = 0; ii < p->Nk; ii++) {
  for (int jj = 0; jj < p->Nc; jj++) {
    vArray->at(p->Ik+ii)[p->Ic+jj] = Vkc;
    vArray->at(p->Ic+jj)[p->Ik+ii] = Vkc;
  }
      }
    }
  }

#ifdef DEBUG
  std::cout << "\nCoupling matrix:\n";
  for (int ii = 0; ii < p->NEQ; ii++) {
    for (int jj = 0; jj < p->NEQ; jj++)
      std::cout << std::scientific << vArray->at(ii)[jj] << " ";
    std::cout << std::endl;
  }
#endif
}

/* Get band index based on flag */
int bandStartIdx(int bandFlag, PARAMETERS * p) {
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
int bandEndIdx(int bandFlag, PARAMETERS * p) {
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
int bandNumStates(int bandFlag, PARAMETERS * p) {
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

void buildParabolicBand(realtype * energies, int n, double bandEdge, int flag, PARAMETERS * p) {
  int s;  // sign +/-
  double m; // mass of electron/hole

  // determine conduction vs. valence band, electron/hole masses
  if (flag == CONDUCTION) {
    s = 1;
    m = p->me;
  }
  else if (flag == VALENCE) {
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
void updateHamiltonian(PARAMETERS * p, realtype t) {
  // TODO unpack NEQ from p

  // get pointer to H
  realtype * H = &(p->H)[0];

  //// first handle torsion
  double torsionValue = 0.0;
  if (p->torsion) {
    // regular (sinusoidal) coupling function
    if (p->torsionSin2) {
      torsionValue = sin2(p->torsionSin2V0, p->torsionSin2V1, p->torsionSin2omega, p->torsionSin2phi, t);
    }
    else {
      torsionValue = p->torsionV.value(t);
    }
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

  // make sparse version of Hamiltonian
  int job [6] = {0, 0, 0, 2, p->NEQ2, 1};
  int info = 0;

  mkl_ddnscsr(&job[0], &(p->NEQ), &(p->NEQ), &(p->H)[0], &(p->NEQ), &(p->H_sp)[0],
      &(p->H_cols)[0], &(p->H_rowind)[0], &info);

  return;
}

#include "params.hpp"

void initialize(PARAMETERS * p) {
  // This function performs error checking on various parameters

  // check torsion parameters, set up torsion spline
  if (p->torsion) {
#ifdef DEBUG
    std::cout << "Torsion is on." << std::endl;
#endif

    // error checking
    if (p->torsionSite > p->Nb) {
      std::cerr << "ERROR: torsion site (" << p->torsionSite
        << ") is larger than number of bridge sites (" << p->Nb << ")." << std::endl;
      exit(-1);
    }
    else if (p->torsionSite < 0) {
      std::cerr << "ERROR: torsion site is less than zero." << std::endl;
      exit(-1);
    }

    if (!p->torsionSin2) {
      if (!fileExists(p->torsionFile)) {
      std::cerr << "ERROR: torsion file " << p->torsionFile << " does not exist." << std::endl;
      exit(-1);
      }

      // create spline
      p->torsionV.readFile(p->torsionFile.c_str());
      if (p->torsionV.getFirstX() != 0.0) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should start at 0.0." << std::endl;
        exit(-1);
      }

      if (p->torsionV.getLastX() < p->tout) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should be >= tout." << std::endl;
        exit(-1);
      }
    }
  }
}

// This function builds up the Hamiltonian, as well as the constituent site
// energy array and coupling array.
void initHamiltonian(PARAMETERS * p) {
  // first read in energies and couplings from files
  p->Nc = numberOfValuesInFile(p->cEnergiesInput.c_str());
  p->Nb = numberOfValuesInFile(p->bEnergiesInput.c_str());

  std::vector<realtype> k_energies (p->Nk);
  std::vector<realtype> c_energies (p->Nc);
  std::vector<realtype> b_energies (p->Nb);
  std::vector<realtype> l_energies (p->Nl);

  // c energies are defined in file
  readVectorFromFile(c_energies, p->cEnergiesInput.c_str(), p->Nc);

  // bridge-dependent parameters
  if (p->bridge_on) {
    if (p->Nb < 1) {
      std::cerr << "\nERROR: bridge_on but no bridge states.  The file b_energies.in is probably empty.\n";
      _exit(-1);
    }

    readVectorFromFile(b_energies, p->bEnergiesInput.c_str(), p->Nb);

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

#ifdef DEBUG
  std::cout << "\nDone reading things from inputs.\n";
#endif

  // assign bulk conduction and valence band energies
  // for RTA, bulk and valence bands have parabolic energies
  if (p->rta) {
    buildParabolicBand(&k_energies[0], p->Nk, p->kBandEdge, CONDUCTION, p);
    buildParabolicBand(&l_energies[0], p->Nl, p->lBandTop, VALENCE, p);
  }
  else {
    buildContinuum(&k_energies[0], p->Nk, p->kBandEdge, p->kBandTop);
    buildContinuum(&l_energies[0], p->Nl, p->kBandEdge - p->valenceBand - p->bulk_gap, p->kBandEdge - p->bulk_gap);
  }
  // calculate band width
  p->kBandWidth = k_energies[p->Nk - 1] - k_energies[0];

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
    p->energies[p->Ik + ii] = k_energies[ii];
  }
  for (int ii = 0; ii < p->Nc; ii++) {
    p->energies[p->Ic + ii] = c_energies[ii];
  }
  for (int ii = 0; ii < p->Nb; ii++) {
    p->energies[p->Ib + ii] = b_energies[ii];
  }
  for (int ii = 0; ii < p->Nl; ii++) {
    p->energies[p->Il + ii] = l_energies[ii];
  }

#ifdef DEBUG
  for (int ii = 0; ii < p->NEQ; ii++) {
    std::cout << "energies[" << ii << "] is " << p->energies[ii] << "\n";
  }
#endif
}

// This function builds up the initial wavefunction coefficients based on inputs
void initWavefunction(PARAMETERS * p) {
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
    if (p->rta) {
      k_coeffs.assign(k_coeffs.size(), 1e-1); // FIXME
    }
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
  p->startWfn.resize(2*p->NEQ, 0.0);

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
#pragma omp parallel for
    for (int ii = 0; ii < p->NEQ; ii++) {
      // diagonal part
      p->startDM[p->NEQ*ii + ii] = pow(p->startWfn[ii],2) + pow(p->startWfn[ii + p->NEQ],2);
      if (p->coherent) {
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
    // No normalization if RTA is on.
    if (!p->rta) {
      summ = 0.0;
      for (int ii = 0; ii < p->NEQ; ii++) {
        // assume here that diagonal elements are all real
        summ += p->startDM[p->NEQ*ii + ii];
      }
      if ( summ == 0.0 ) {
        std::cerr << "\nFATAL ERROR [populations]: total population is 0!\n";
        _exit(-1);
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
    }

    // Error checking for total population; recount population first
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += p->startDM[p->NEQ*ii + ii];
    }
    if ( fabs(summ-1.0) > 1e-12  && (!p->rta)) {
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



int main (int argc, char * argv[]) {



  //// DECLARING VARIABLES


  // Struct of parameters
  PARAMETERS p;
  // CVode variables
  void * cvode_mem = NULL;                      // pointer to block of CVode memory
  N_Vector y, yout;                     // arrays of populations

  // arrays for energetic parameters
  std::vector< std::vector<realtype> > V;

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
  std::map<const std::string, bool> outs;       // map of output file names to bool

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

  //// ASSIGN PARAMETERS FROM INPUT FILE


  // ---- TODO create output directory if it does not exist ---- //
  flag = mkdir(p.outputDir.c_str(), 0755);

  assignParams(p.inputFile.c_str(), &p);

  // Decide which output files to make
#ifdef DEBUG
  std::cout << "Assigning outputs as specified in " << p.inputFile << "\n";
#endif
  assignOutputs(p.inputFile.c_str(), outs, &p);
  p.outs = outs;

#ifdef DEBUG
  // print out which outputs will be made
  for (std::map<const std::string, bool>::iterator it = outs.begin(); it != outs.end(); it++) {
    std::cout << "Output file: " << it->first << " will be created.\n";
  }
#endif

  // OPEN LOG FILE; PUT IN START TIME //
  if (isOutput(outs, "log.out")) {
    log = fopen("log.out", "w");                        // note that this file is closed at the end of the program
  }
  time(&startRun);
  currentTime = localtime(&startRun);
  if (isOutput(outs, "log.out")) {
    fprintf(log, "Run started at %s\n", asctime(currentTime));
  }

  if (isOutput(outs, "log.out")) {
    // make a note about the laser intensity.
    fprintf(log,"The laser intensity is %.5e W/cm^2.\n\n",pow(p.pumpAmpl,2)*3.5094452e16);
  }


  //// READ DATA FROM INPUTS
  // TODO replace with initialize function
  initHamiltonian(&p);
  initWavefunction(&p);

  //// PREPROCESS DATA FROM INPUTS


  // set number of processors for OpenMP
  omp_set_num_threads(p.nproc);
  mkl_set_num_threads(p.nproc);
#ifdef DEBUG
  std::cout << "\nTotal number of states: " << p.NEQ << std::endl;
  std::cout << p.Nk << " bulk, " << p.Nc << " QD, " << p.Nb << " bridge, " << p.Nl << " bulk VB.\n";
#endif
  // assign times.
  p.times.resize(p.numOutputSteps+1);
  for (int ii = 0; ii <= p.numOutputSteps; ii++) {
    p.times[ii] = float(ii)/p.numOutputSteps*p.tout;
  }


  //// ASSIGN COUPLING CONSTANTS

  buildCoupling(&V, &p, outs);

  if (isOutput(outs, "log.out")) {
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

  //// BUILD HAMILTONIAN


  // //TODO TODO
#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif
  realtype * H = NULL;
  H = new realtype [p.NEQ2];
  for (int ii = 0; ii < p.NEQ2; ii++) {
    H[ii] = 0.0;
  }
  buildHamiltonian(H, p.energies, &V, &p);
  // add Hamiltonian to p
  p.H.resize(p.NEQ2);
  for (int ii = 0; ii < p.NEQ2; ii++) {
    p.H[ii] = H[ii];
  }
  // create sparse version of H
  p.H_sp.resize(p.NEQ2);
  p.H_cols.resize(p.NEQ2);
  p.H_rowind.resize(p.NEQ2 + 1);
  int job [6] = {0, 0, 0, 2, p.NEQ2, 1};
  int info = 0;

  mkl_ddnscsr(&job[0], &(p.NEQ), &(p.NEQ), &(p.H)[0], &(p.NEQ), &(p.H_sp)[0],
      &(p.H_cols)[0], &(p.H_rowind)[0], &info);


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
  makePlots(outs, &p);

  // only do propagation if not just making plots
  if (! p.justPlots) {
    // Make outputs independent of time propagation
    computeGeneralOutputs(outs, &p);
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


    //// MAKE FINAL OUTPUTS


    // finalize log file //
    time(&endRun);
    currentTime = localtime(&endRun);
    if (isOutput(outs, "log.out")) {
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
      computeWfnOutput(wfnt, outs, &p);
    }
    else {
      computeDMOutput(dmt, outs, &p);
    }
#ifdef DEBUG
    std::cout << "done computing outputs" << std::endl;
#endif

    // do analytical propagation
    if (p.analytical && (! p.bridge_on)) {
      computeAnalyticOutputs(outs, &p);
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
  delete [] H;
  if (p.wavefunction) {
    delete [] wfnt;
  }
  else {
    delete [] dmt;
  }

  std::cout << "whoo" << std::endl;

  return 0;
}
