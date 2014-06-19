#include <params.hpp>

//#define DEBUG_BUILDHAMILTONIAN

/* builds a Hamiltonian from site energies and couplings. */
void PARAMETERS::buildHamiltonian() {
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