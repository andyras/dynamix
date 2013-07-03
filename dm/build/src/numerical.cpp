#include "numerical.h"

/* returns the number of numbers in a file.  This way, it doesn't matter if
 * they are one per line or multiple per line.
 */
int numberOfValuesInFile(const char * nameOfFile) {
 FILE * input;
 double value;
 int numberOfValues = 0;

 input = fopen(nameOfFile, "r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &value) != EOF) { numberOfValues++; }
  if (numberOfValues == 0 ) {
   fprintf(stderr, "WARNING: file %s is empty.\n", nameOfFile);
  }
 }
 else {
  fprintf(stderr, "WARNING [numberOfValuesInFile]: file %s does not exist.\n", nameOfFile);
  return -1;
 }

 return numberOfValues;
}

/* reads in the values from file; returns an array the length of the number of 
 * numbers in the file
 */
void readArrayFromFile(realtype * array, const char * nameOfFile, int numberOfValues) {
 FILE * input;
 int i = 0;

 input = fopen(nameOfFile,"r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &array[i]) != EOF && i < numberOfValues) {
   i++;
  }
 }
 else {
  fprintf(stderr, "ERROR [readArrayFromFile]: file %s does not exist.\n", nameOfFile);
 }

 fclose(input);
}

/* Returns an array of length n with all values set to initializeValue. */
void initializeArray(realtype * array, int n, realtype initializeValue) {
#ifdef DEBUG
 cout << "initializeValue is " << initializeValue << endl;
#endif

 int i;

 for (i = 0; i < n; i++) {
  array[i] = initializeValue;
 }
}

/* builds energies for a quasicontinuum (evenly spaced) */
void buildContinuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
 
 int i;

 Energies[0] = BandEdge;	// the bottom of the conduction band is set
 
 // loop over the remaining states.  This way the top of the band will be at BandTop
 for (i = 1; i < numberOfStates; i++) {
  Energies[i] = Energies[i-1] + (BandTop-BandEdge)/(numberOfStates-1);
 }
}

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void buildKPops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk) {

 int i;

 for (i = 0; i < Nk; i++) {
  kPops[i] = sqrt(1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.01)*3.185e5/(temp))));
#ifdef DEBUG
 cout << "\nk population at state " << i << " is: "
      << sqrt(1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.01)*3.185e5/(temp))));
#endif
 }
#ifdef DEBUG
 cout << endl;
#endif
}

/* populates a set of states according to a Gaussian distribution. */
void buildKPopsGaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk) {

 int i;

 for (i = 0; i < Nk; i++) {
  kPops[i] = sqrt((1/(sigma*sqrt(2*3.1415926535)))*exp(-pow((kEnergies[i]-(kBandEdge+mu)),2)/(2*pow(sigma,2))));
#ifdef DEBUG
 cout << "\nk population at state " << i << " is: "
      << sqrt((1/(sigma*sqrt(2*3.1415926535)))*exp(-pow((kEnergies[i]-(kBandEdge+mu)),2)/(2*pow(sigma,2))));
#endif
 }
#ifdef DEBUG
 cout << endl;
#endif
}

/* returns the coupling as a function of energy E given that the middle of the
 * band is at position mid.
 * Eq. 31 in Ramakrishna et al., JCP 2001, 115, 2743-2756
 */
double parabolicV(double Vee, double E, double bandEdge, double bandTop) {
// set band edge as zero energy
E = E - bandEdge;
double mid = (bandTop - bandEdge)/2.0;

#ifdef DEBUG
 fprintf(stdout, "coupling at (E - band edge) = %.9e: %.9e\n", E, Vee*sqrt(sqrt(pow(mid,2) - pow((E - mid),2))/mid));
#endif
 // test whether at the very top or bottom of band
 if (abs(abs(E - mid) - mid) < 1e-10) {
#ifdef DEBUG
  fprintf(stdout, "(at bottom or top of band edge, returning 0.0)\n");
#endif
  return 0.0;
 }
 return Vee*sqrt(sqrt(pow(mid,2) - pow((E - mid),2))/mid);
}

/* gives the value of a laser pulse (electric field) at time t */
realtype pump(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase) {
 double sigma = pumpFWHM/2.35482005;
 return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total) {
 int i;
 realtype summ = 0;

 for (i = 0; i < NV_LENGTH_S(nv); i++) {
  summ += (NV_Ith_S(nv, i)*NV_Ith_S(nv, i));
 }
 summ = sqrt(summ);
 for (i = 0; i < NV_LENGTH_S(nv); i++) {
  NV_Ith_S(nv, i) = total*NV_Ith_S(nv,i)/summ;
 }

 return 0;
}

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep) {
 int i;		// counter

 if (inputLength < 6 ) {
  fprintf(stderr, "ERROR [Derivative]: array has length less than 6 elements, cannot proceed");
  return -1;
 }

 for (i = 2; i < inputLength-3; i++) {
  outputArray[i-2] = (2* inputArray[i+3]
		     -15*inputArray[i+2]
		     +60*inputArray[i+1]
		     -20*inputArray[i]
		     -30*inputArray[i-1]
		     +3 *inputArray[i-2])/(60*timestep);
 }

 return 0;
}

/* Riemann sum of an array (values) at time points (time).
 * Does not assume equal spacing in time.
 */
realtype integrateArray(realtype * values, realtype * time, int num) {
 int i;
 realtype riemann = 0;

 for (i = 0; i < num-1; i++)
  riemann += (values[i+1] + values[i])*(time[i+1]-time[i])/2;

 return riemann;
}

/* Returns maximum element in an array. */
realtype findArrayMaximum(realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];

 for (i = 1; i < num; i++)
  if (inputArray[i] > currentMax)
   currentMax = inputArray[i];

 return currentMax;
}

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype findFirstArrayMaximum(realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax)
   currentMax = inputArray[i];
  if (inputArray[i] < currentMax)
   break;
 }

 return currentMax;
}

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int findFirstArrayMaximumIndex(realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];
 int currentMax_index = 0;

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax) {
   currentMax = inputArray[i];
   currentMax_index = i;
  }
  if (inputArray[i] < currentMax)
   break;
 }

 return currentMax_index;
}

/* returns index of first maximum in an array. */
int findArrayMaximumIndex(realtype * inputArray, int num) {
 int i;
 realtype currentMax = inputArray[0];
 int currentMax_index = 0;

 for (i = 1; i < num; i++) {
  if (inputArray[i] > currentMax) {
   currentMax = inputArray[i];
   currentMax_index = i;
  }
 }

 return currentMax_index;
}

/* assign coupling constants to global array V */
void buildCoupling (realtype ** vArray, struct PARAMETERS * p,
                    std::map<std::string, bool> &outs) {
 
 int i, j;	// counters
 double Vkc;	// coupling between bulk and QD
 double Vkb1;	// coupling between bulk and first bridge
 double VbNc;	// coupling between last bridge and QD

 // initialize the coupling array
 for (i = 0; i < p->NEQ; i++) {
  for (j = 0; j < p->NEQ; j++) {
   vArray[i][j] = 0.0;
  }
 }

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
   for (i = 0; i < p->Nk; i++) {
    vArray[p->Ik+i][p->Ib] = parabolicV(Vkb1, p->energies[p->Ik+i], p->kBandEdge, p->kBandTop);
    vArray[p->Ib][p->Ik+i] = parabolicV(Vkb1, p->energies[p->Ik+i], p->kBandEdge, p->kBandTop);
   }
  }
  else {
   for (i = 0; i < p->Nk; i++) {
    vArray[p->Ik+i][p->Ib] = Vkb1;
    vArray[p->Ib][p->Ik+i] = Vkb1;
   }
  }
   
  // coupling between bN and c
  if ((p->scale_brqd) && (p->Nc > 1)) {
   VbNc = p->Vbridge[p->Nb]/sqrt(p->Nc-1);
  }
  else {
   VbNc = p->Vbridge[p->Nb];
  }
  for (i = 0; i < p->Nc; i++) {
   vArray[p->Ic+i][p->Ib+p->Nb-1] = VbNc;
   vArray[p->Ib+p->Nb-1][p->Ic+i] = VbNc;
  }
  
  // coupling between bridge states
  for (i = 0; i < p->Nb - 1; i++) {
   vArray[p->Ib+i][p->Ib+i+1] = p->Vbridge[i+1];
   vArray[p->Ib+i+1][p->Ib+i] = p->Vbridge[i+1];
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
   for (i = 0; i < p->Nk; i++) {
    for (j = 0; j < p->Nc; j++) {
     vArray[p->Ik+i][p->Ic+j] = parabolicV(Vkc, p->energies[p->Ik+i], p->kBandEdge, p->kBandTop);
     vArray[p->Ic+j][p->Ik+i] = parabolicV(Vkc, p->energies[p->Ik+i], p->kBandEdge, p->kBandTop);
    }
   }
  }
  else {
   for (i = 0; i < p->Nk; i++) {
    for (j = 0; j < p->Nc; j++) {
     vArray[p->Ik+i][p->Ic+j] = Vkc;
     vArray[p->Ic+j][p->Ik+i] = Vkc;
    }
   }
  }
 }

#ifdef DEBUG
 cout << "\nCoupling matrix:\n";
 for (i = 0; i < p->NEQ; i++) {
  for (j = 0; j < p->NEQ; j++)
   cout << scientific << vArray[i][j] << " ";
  cout << endl;
 }
#endif

 FILE * couplings;
 if (outs["couplings.out"]) {
  couplings = fopen("couplings.out","w");
  for (i = 0; i < p->NEQ; i++) {
   for (j = 0; j < p->NEQ; j++) {
    fprintf(couplings,"%.7g ",vArray[i][j]);
   }
   fprintf(couplings,"\n");
  }
  fclose(couplings);
 }
}

/* builds a Hamiltonian from site energies and couplings. */
void buildHamiltonian(realtype * H, realtype * energy, realtype ** V, struct PARAMETERS * p) {
 // indices
 int idx1, idx2;
 int N = p->NEQ;
 
#ifdef DEBUG
  fprintf(stderr, "Assigning diagonal elements of Hamiltonian.\n");
#endif
  for (int ii = 0; ii < N; ii++) {
   // diagonal
   H[ii*N + ii] = energy[ii];
#ifdef DEBUG
   cout << "diagonal element " << ii << " of H is " << energy[ii] << "\n";
#endif
  }

 if (p->bridge_on) {
  // assign bulk-bridge coupling
#ifdef DEBUG
  fprintf(stderr, "Assigning bulk-bridge coupling elements in Hamiltonian.\n");
#endif
  idx2 = p->Ib;
  for (int ii = 0; ii < p->Nk; ii++) {
   idx1 = p->Ik + ii;
#ifdef DEBUG2
   fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
   fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
   fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
   H[idx1*N + idx2] = V[idx1][idx2];
   H[idx2*N + idx1] = V[idx2][idx1];
  }
  // assign bridge-bridge couplings
#ifdef DEBUG
  fprintf(stderr, "Assigning bridge-bridge coupling elements in Hamiltonian.\n");
#endif
  for (int ii = 1; ii < p->Nb; ii++) {
   idx1 = p->Ib + ii;
   idx2 = p->Ib+ ii + 1;
#ifdef DEBUG2
   fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
   fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
   fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
   H[idx1*N + idx2] = V[idx1][idx2];
   H[idx2*N + idx1] = V[idx2][idx1];
  }
  // assign bridge-QD coupling
#ifdef DEBUG
  fprintf(stderr, "Assigning bridge-QD coupling elements in Hamiltonian.\n");
#endif
  idx2 = p->Ib + p->Nb - 1;
  for (int ii = 0; ii < p->Nc; ii++) {
   idx1 = p->Ic + ii;
#ifdef DEBUG2
   fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
   fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
   fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
   H[idx1*N + idx2] = V[idx1][idx2];
   H[idx2*N + idx1] = V[idx2][idx1];
  }
 }
 // no bridge
 else {
  // assign bulk-QD coupling
#ifdef DEBUG
  fprintf(stderr, "Assigning bulk-QD coupling elements in Hamiltonian.\n");
#endif
  for (int ii = 0; ii < p->Nk; ii++) {
   idx1 = p->Ik + ii;
   for (int jj = 0; jj < p->Nc; jj++) {
    idx2 = p->Ic + jj;
#ifdef DEBUG2
    fprintf(stderr, "H[%d*%d + %d] = ", idx1, N, idx2);
    fprintf(stderr, "V[%d][%d] = ", idx1, idx2);
    fprintf(stderr, "%e\n", V[idx1][idx2]);
#endif
    H[idx1*N + idx2] = V[idx1][idx2];
    H[idx2*N + idx1] = V[idx2][idx1];
   }
  }
 }
}

/* Updates \rho(t) at each time step. */
void updateDM(N_Vector dm, realtype * dmt, int timeStep, struct PARAMETERS * p) {
 for (int ii = 0; ii < p->NEQ2; ii++) {
  dmt[2*p->NEQ2*timeStep + ii] = NV_Ith_S(dm, ii);
  dmt[2*p->NEQ2*timeStep + ii + p->NEQ2] = NV_Ith_S(dm, ii + p->NEQ2);
 }

 return;
}

