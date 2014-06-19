#include <params.hpp>

#define DEBUG_BUILDCOUPLING
//#define DEBUG_BUILDHAMILTONIAN

/* assign coupling constants to global array V */
void Params::buildCoupling () {

  double Vkc; // coupling between bulk and QD
  double Vkb1;  // coupling between bulk and first bridge
  double VbNc;  // coupling between last bridge and QD

  // initialize the coupling array
  V.resize(NEQ);
  for (int ii = 0; ii < V.size(); ii++) {
    V.at(ii).assign(NEQ, 0.0);
  }

std::cout << "\n\n\nWHOOOOOT\n\n\n";

#ifdef DEBUG_BUILDCOUPLING
  if (bridge_on) {
    for (int ii = 0; ii < Nb + 1; ii++) {
      std::cout << "Vbridge[" << ii << "] is ";
      std::cout << Vbridge[ii] << "\n";
    }
  }
#endif

  // bridge
  if (bridge_on) {
    // coupling between k and b1
    if ((scale_bubr) && (Nk > 1)) {
      Vkb1 = sqrt(Vbridge[0]*(kBandTop-kBandEdge)/(Nk-1));
    }
    else {
      Vkb1 = Vbridge[0];
    }
    if (parabolicCoupling) {
      for (int ii = 0; ii < Nk; ii++) {
        V.at(Ik+ii)[Ib] = parabolicV(Vkb1, energies[Ik+ii], kBandEdge, kBandTop);
        V.at(Ib)[Ik+ii] = parabolicV(Vkb1, energies[Ik+ii], kBandEdge, kBandTop);
      }
    }
    else {
      for (int ii = 0; ii < Nk; ii++) {
        V.at(Ik+ii)[Ib] = Vkb1;
        V.at(Ib)[Ik+ii] = Vkb1;
      }
    }

    // coupling between bN and c
    if ((scale_brqd) && (Nc > 1)) {
      VbNc = Vbridge[Nb]/sqrt(Nc-1);
    }
    else {
      VbNc = Vbridge[Nb];
    }
    for (int ii = 0; ii < Nc; ii++) {
      V.at(Ic+ii)[Ib+Nb-1] = VbNc;
      V.at(Ib+Nb-1)[Ic+ii] = VbNc;
    }

    // coupling between bridge states
    for (int ii = 0; ii < Nb - 1; ii++) {
      V.at(Ib+ii)[Ib+ii+1] = Vbridge[ii+1];
      V.at(Ib+ii+1)[Ib+ii] = Vbridge[ii+1];
    }
  }
  // no bridge
  else {
    // scaling
    if ((scale_buqd) && (Nk > 1)) {
      Vkc = sqrt(Vnobridge[0]*(kBandTop-kBandEdge)/(Nk-1));
    }
    else {
      Vkc = Vnobridge[0];
    }

    // parabolic coupling of bulk band to QD
    if (parabolicCoupling) {
      for (int ii = 0; ii < Nk; ii++) {
        for (int jj = 0; jj < Nc; jj++) {
          V.at(Ik+ii)[Ic+jj] = parabolicV(Vkc, energies[Ik+ii], kBandEdge, kBandTop);
          V.at(Ic+jj)[Ik+ii] = parabolicV(Vkc, energies[Ik+ii], kBandEdge, kBandTop);
        }
      }
    }
    else {
      for (int ii = 0; ii < Nk; ii++) {
        for (int jj = 0; jj < Nc; jj++) {
          V.at(Ik+ii)[Ic+jj] = Vkc;
          V.at(Ic+jj)[Ik+ii] = Vkc;
        }
      }
    }
  }

#ifdef DEBUG
  std::cout << "\nCoupling matrix:\n";
  for (int ii = 0; ii < NEQ; ii++) {
    for (int jj = 0; jj < NEQ; jj++)
      std::cout << std::scientific << V.at(ii)[jj] << " ";
    std::cout << std::endl;
  }
#endif
}

/* builds a Hamiltonian from site energies and couplings. */
void Params::buildHamiltonian() {
  // indices
  int idx1, idx2;
  int N = NEQ;

  H.assign(N*N, 0.0);

#ifdef DEBUG_BUILDHAMILTONIAN
  fprintf(stderr, "Assigning diagonal elements of Hamiltonian.\n");
#endif
  for (int ii = 0; ii < N; ii++) {
    // diagonal
    H[ii*N + ii] = energies[ii];
#ifdef DEBUG_BUILDHAMILTONIAN
    std::cout << "diagonal element " << ii << " of H is " << energies[ii] << "\n";
#endif
  }

  if (bridge_on) {
    // assign bulk-bridge coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-bridge coupling elements in Hamiltonian.\n");
#endif
    idx2 = Ib;
    for (int ii = 0; ii < Nk; ii++) {
      idx1 = Ik + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V.at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V.at(idx1)[idx2];
      H[idx2*N + idx1] = V.at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bulk-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-bridge couplings
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-bridge coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < (Nb-1); ii++) {
      idx1 = Ib + ii;
      idx2 = Ib+ ii + 1;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V.at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V.at(idx1)[idx2];
      H[idx2*N + idx1] = V.at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bridge-bridge coupling elements in Hamiltonian.\n");
    // assign bridge-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bridge-QD coupling elements in Hamiltonian.\n");
#endif
    idx2 = Ib + Nb - 1;
    for (int ii = 0; ii < Nc; ii++) {
      idx1 = Ic + ii;
#ifdef DEBUG_BUILDHAMILTONIAN
      fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
      fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
      fprintf(stderr, "%e\n", V.at(idx1)[idx2]);
#endif
      H[idx1*N + idx2] = V.at(idx1)[idx2];
      H[idx2*N + idx1] = V.at(idx2)[idx1];
    }
    fprintf(stderr, "Done assigning bridge-QD coupling elements in Hamiltonian.\n");
  }
  // no bridge
  else {
    // assign bulk-QD coupling
#ifdef DEBUG_BUILDHAMILTONIAN
    fprintf(stderr, "Assigning bulk-QD coupling elements in Hamiltonian.\n");
#endif
    for (int ii = 0; ii < Nk; ii++) {
      idx1 = Ik + ii;
      for (int jj = 0; jj < Nc; jj++) {
  idx2 = Ic + jj;
#ifdef DEBUG_BUILDHAMILTONIAN
  fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
  fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
  fprintf(stderr, "%e\n", V.at(idx1)[idx2]);
#endif
  H[idx1*N + idx2] = V.at(idx1)[idx2];
  H[idx2*N + idx1] = V.at(idx2)[idx1];
      }
    }
  }

#ifdef DEBUG
  std::cout << "\nTotal number of states: " << NEQ << std::endl;
  std::cout << Nk << " bulk, " << Nc << " QD, " << Nb << " bridge, " << Nl << " bulk VB.\n";
#endif
  // assign times.
  times.resize(numOutputSteps+1);
  for (int ii = 0; ii <= numOutputSteps; ii++) {
    times[ii] = float(ii)/numOutputSteps*tout;
  }
}