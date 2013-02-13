#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <time.h>
#include <numeric>
#include <complex>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>

/* DEBUG compiler flag: turn on to generate basic debug outputs.         */
// #define DEBUG
// DEBUG2 flag: turn on for more numerical output
//#define DEBUG2
/* DANGER! Only turn on DEBUGf for small test runs, otherwise output is       */
/* enormous (many GB).  This flag turns on debug output within the f          */
/* function.                                                                  */
//#define DEBUGf
/* This flag is debuggery related to checking against Sai's code.             */
//#define DEBUG_SAI

using namespace std;

// structure to hold complex numbers for linear algebra routines
typedef struct {
 double re;
 double im;
} complex16;

// GLOBAL VARIABLES GO HERE //
void * cvode_mem;			// pointer to block of CVode memory
realtype * user_data;
N_Vector y, yout;			// arrays of populations
int Nk;				// number of each type of state
int Nc;
int Nb;
int Nl;
int N_vib;				// number of vibronic states
int Ik;				// index starters for each type of state
int Ic;
int Ib;
int Il;
int Ik_vib;				// index starters for each type of vibronic state
int Ic_vib;
int Ib_vib;
int Il_vib;
int NEQ;				// total number of states/equations
int NEQ_vib;
int numOutputSteps;			// number of timesteps
realtype k_bandedge;			// lower edge of bulk conduction band
realtype k_bandtop;			// upper edge of bulk conduction band
double E_vib;				// vibrational energy
double gkc;				// g factor between k and c states
double gkb;				// g factor between k and b states
double gbc;				// g factor between b and c states
double gbb;				// g factor between b states
double muLK;                           // transition dipole moment from l to k (energy a.u.)
double pumpFWHM;                       // FWHM of pump pulse (time a.u.)
double pumpPeak;                       // time of peak of pump pulse (a.u.)
double pumpFreq;                       // frequency of pump pulse (energy a.u.)
double pumpAmpl;                       // intensity of pump pulse (electric field a.u.)
double pumpPhase;                      // pump pulse phase (in units of radians)
realtype ** V;				// pointer to k-c coupling constants
realtype ** FCkc;			// Franck-Condon factors
realtype ** FCkb;
realtype ** FCbb;
realtype ** FCbc;
realtype * energy;
realtype * Vbridge;			// pointer to array of bridge coupling constants.
					// first element [0] is Vkb1, last [Nb] is VcbN
realtype * Vnobridge;			// coupling constant when there is no bridge
bool bulk_FDD = 0;			// switches for starting conditions
bool bulk_Gauss = 0;
bool bulk_constant = 0;
bool qd_pops = 0;
bool laser_on = 0;
bool parabolicCoupling = 0;
bool scale_bubr = 0;
bool scale_brqd = 0;
bool scale_buqd = 0;
bool scale_laser = 0;
bool bridge_on = 0;
bool random_phase = 0;
int random_seed = 0;
// END GLOBAL VARIABLES
#ifdef DEBUG_SAI
double last_t = -1.0;				// keeps track of last time for which debuggery was printed
#endif

int Number_of_values (const char * nameOfFile) {
// returns the number of numbers in a file.  This way, it doesn't matter if
// they are one per line or multiple per line.
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
  fprintf(stderr, "WARNING [Number_of_values]: file %s does not exist.\n", nameOfFile);
  return -1;
 }

 return numberOfValues;
}

void Read_array_from_file (realtype * array, const char * nameOfFile, int numberOfValues) {
// reads in the values from file; returns an array the length of the number of 
// numbers in the file

 FILE * input;
 int i = 0;

 input = fopen(nameOfFile,"r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &array[i]) != EOF && i < numberOfValues) {
   i++;
  }
 }
 else {
  fprintf(stderr, "ERROR [Read_array_from_file]: file %s does not exist.\n", nameOfFile);
 }

 fclose(input);
}

void Initialize_array(realtype * array, int n, realtype initializeValue) {
// Returns an array of length n with all values set to initializeValue. //
#ifdef DEBUG
 cout << "initializeValue is " << initializeValue << endl;
#endif

 int i;

 for (i = 0; i < n; i++) {
  array[i] = initializeValue;
 }
}


void Build_continuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
// builds energies for a quasicontinuum (evenly spaced)
 
 int i;

 Energies[0] = BandEdge;	// the bottom of the conduction band is set
 
 // loop over the remaining states.  This way the top of the band will be at BandTop
 for (i = 1; i < numberOfStates; i++) {
  Energies[i] = Energies[i-1] + (BandTop-BandEdge)/(numberOfStates-1);
 }
}


void Build_k_pops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp) {
// populates a set of states according to a Fermi-Dirac distribution.
// I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
// below the lowest-energy state in the set of states being populated.

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


void Build_k_pops_Gaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu) {
// populates a set of states according to a Gaussian distribution.

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


void Build_Franck_Condon_factors (realtype ** FCmat, double g, int numM, int numN) {
// creates a matrix of Franck-Condon factors between two displaced harmonic
// oscillators.

 int m, n;

 FCmat[0][0] = exp(-pow(g,2)/2);	// first element
 for (m = 1; m < numM; m++)		// first column
  FCmat[m][0] = g/sqrt((double)m)*FCmat[m-1][0];
 for (n = 1; n < numN; n++)		// first row
  FCmat[0][n] = -1*g/sqrt((double)n)*FCmat[0][n-1];
 for (m = 1; m < numM; m++)		// recursion for the rest of the matrix
  for (n = 1; n < numN; n++)
   FCmat[m][n] = g/sqrt((double)m)*FCmat[m-1][n] + sqrt((double)n/(double)m)*FCmat[m-1][n-1];

}

double parabolicV(double Vee, double E, double bandEdge, double bandTop) {
// returns the coupling as a function of energy E given that the middle of the
// band is at position mid.
// Eq. 31 in Ramakrishna et al., JCP 2001, 115, 2743-2756

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

void buildCoupling (realtype ** vArray, int dim, realtype kBandEdge, realtype kBandTop, realtype * energy) {
// assign coupling constants to global array V
 
 int i, j;	// counters
 double Vkc;	// coupling between bulk and QD
 double Vkb1;	// coupling between bulk and first bridge
 double VbNc;	// coupling between last bridge and QD

 // initialize the coupling array
 for (i = 0; i < dim; i++) {
  for (j = 0; j < dim; j++) {
   vArray[i][j] = 0.0;
  }
 }

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
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = parabolicV(Vkb1, energy[Ik+i], kBandEdge, kBandTop);
    vArray[Ib][Ik+i] = parabolicV(Vkb1, energy[Ik+i], kBandEdge, kBandTop);
   }
  }
  else {
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = Vkb1;
    vArray[Ib][Ik+i] = Vkb1;
   }
  }
   
  // coupling between bN and c
  if ((scale_brqd) && (Nc > 1)) {
   VbNc = Vbridge[Nb]/sqrt(Nc-1);
  }
  else {
   VbNc = Vbridge[Nb];
  }
  for (i = 0; i < Nc; i++) {
   vArray[Ic+i][Ib+Nb-1] = VbNc;
   vArray[Ib+Nb-1][Ic+i] = VbNc;
  }
  
  // coupling between bridge states
  for (i = 0; i < Nb - 1; i++) {
   vArray[Ib+i][Ib+i+1] = Vbridge[i+1];
   vArray[Ib+i+1][Ib+i] = Vbridge[i+1];
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
#ifdef DEBUG_SAI
  Vkc = Vnobridge[0]/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);
#endif

  // parabolic coupling of bulk band to QD
  if (parabolicCoupling) {
   for (i = 0; i < Nk; i++) {
    for (j = 0; j < Nc; j++) {
     vArray[Ik+i][Ic+j] = parabolicV(Vkc, energy[Ik+i], kBandEdge, kBandTop);
     vArray[Ic+j][Ik+i] = parabolicV(Vkc, energy[Ik+i], kBandEdge, kBandTop);
    }
   }
  }
  else {
   for (i = 0; i < Nk; i++) {
    for (j = 0; j < Nc; j++) {
     vArray[Ik+i][Ic+j] = Vkc;
     vArray[Ic+j][Ik+i] = Vkc;
    }
   }
  }
 }

#ifdef DEBUG
 cout << "\nCoupling matrix:\n";
 for (i = 0; i < dim; i++) {
  for (j = 0; j < dim; j++)
   cout << scientific << vArray[i][j] << " ";
  cout << endl;
 }
#endif

 FILE * couplings;
 couplings = fopen("couplings.out","w");
 for (i = 0; i < dim; i++) {
  for (j = 0; j < dim; j++) {
   fprintf(couplings,"%.7g ",vArray[i][j]);
  }
  fprintf(couplings,"\n");
 }
 fclose(couplings);
 
}

realtype pump(realtype t) {
// gives the value of a laser pulse (electric field) at time t
 double sigma = pumpFWHM/2.35482005;
 return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}


int f(realtype t, N_Vector y, N_Vector ydot, void * data) {
// gives f(y,t) for CVODE
 
 int i, j, n, m;
 int IkRe, IkIm, IcRe, IcIm, IlRe, IlIm;
 realtype sinn;
 realtype coss;
 realtype Vee;
 // initialization
 for (i = 0; i < 2*NEQ_vib; i++)
  NV_Ith_S(ydot, i) = 0;

 if (laser_on) {
  realtype pumpTerm = 0.0;
  if ((scale_laser) && (Nk > 1)) {
   // pumpTerm gives the strength of the pump's interaction, accounting for scaling.
   pumpTerm = muLK*pump(t)*sqrt((k_bandtop - k_bandedge)/(Nk - 1));
  }
  else {
   // pumpTerm gives the strength of the pump's interaction, not accounting for scaling.
   pumpTerm = muLK*pump(t);
  }
  // pump pulse coupling l and k states
  for (i = 0; i < Nk; i++)
   for (j = 0; j < Nl; j++)
    for (n = 0; n < N_vib; n++)
     for (m = 0; m < N_vib; m++) {
      IkRe = Ik_vib + i*N_vib + n;
      IkIm = IkRe + NEQ_vib;
      IlRe = Il_vib + j*N_vib + m;
      IlIm = IlRe + NEQ_vib;
      coss = cos((energy[IkRe] - energy[IlRe])*t);
      sinn = sin((energy[IkRe] - energy[IlRe])*t);
      NV_Ith_S(ydot, IkRe) += pumpTerm*(coss*NV_Ith_S(y, IlIm) + sinn*NV_Ith_S(y, IlRe)); // k Re
      NV_Ith_S(ydot, IkIm) += pumpTerm*(sinn*NV_Ith_S(y, IlIm) - coss*NV_Ith_S(y, IlRe)); // k Im
      NV_Ith_S(ydot, IlRe) += pumpTerm*(coss*NV_Ith_S(y, IkIm) - sinn*NV_Ith_S(y, IkRe)); // l Re
      NV_Ith_S(ydot, IlIm) -= pumpTerm*(sinn*NV_Ith_S(y, IkIm) + coss*NV_Ith_S(y, IkRe)); // l Im
#ifdef DEBUGf
      cout << endl << "IkRe " << IkRe << " IkIm " << IkIm << " IlRe " << IcRe << " IlIm ";
      cout << IcIm << " V " << Vee << " cos " << coss << " sin " << sinn << " t " << t << endl;
#endif
     }
 }

 if (bridge_on) {					// bridge
  int IbRe, IbIm, IBRe, IBIm;
  realtype Vkb = V[Ik][Ib];
  realtype Vbc = V[Ic][Ib+Nb-1];
  for (i = 0; i < Nk; i++)				// Vkb
   for (n = 0; n < N_vib; n++)
    for (m = 0; m < N_vib; m++) {
     IkRe = Ik_vib + i*N_vib + n;
     IkIm = IkRe + NEQ_vib;
     IbRe = Ib_vib + m;
     IbIm = IbRe + NEQ_vib;
     Vee = Vkb*FCkb[n][m];
     coss = Vee*cos((energy[IkRe] - energy[IbRe])*t);
     sinn = Vee*sin((energy[IkRe] - energy[IbRe])*t);
     NV_Ith_S(ydot, IkRe) += (coss*NV_Ith_S(y, IbIm) + sinn*NV_Ith_S(y, IbRe));	// k Re
     NV_Ith_S(ydot, IkIm) += (sinn*NV_Ith_S(y, IbIm) - coss*NV_Ith_S(y, IbRe));	// k Im
     NV_Ith_S(ydot, IbRe) += (coss*NV_Ith_S(y, IkIm) - sinn*NV_Ith_S(y, IkRe));	// b1 Re
     NV_Ith_S(ydot, IbIm) -= (sinn*NV_Ith_S(y, IkIm) + coss*NV_Ith_S(y, IkRe));	// b1 Im
#ifdef DEBUGf
     cout << endl << "IkRe " << IkRe << " IkIm " << IkIm << " IbRe " << IbRe << " IbIm ";
     cout << IbIm << " Vkb " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
#endif
    }
  for (i = 0; i < Nc; i++)				// Vcb
   for (n = 0; n < N_vib; n++)
    for (m = 0; m < N_vib; m++) {
     IcRe = Ic_vib + i*N_vib + n;
     IcIm = IcRe + NEQ_vib;
     IbRe = Ib_vib + (Nb-1)*N_vib + m;
     IbIm = IbRe + NEQ_vib;
     Vee = Vbc*FCbc[m][n];
     coss = Vee*cos((energy[IcRe] - energy[IbRe])*t);
     sinn = Vee*sin((energy[IcRe] - energy[IbRe])*t);
     NV_Ith_S(ydot, IcRe) += (coss*NV_Ith_S(y, IbIm) + sinn*NV_Ith_S(y, IbRe));	// c Re
     NV_Ith_S(ydot, IcIm) += (sinn*NV_Ith_S(y, IbIm) - coss*NV_Ith_S(y, IbRe));	// c Im
     NV_Ith_S(ydot, IbRe) += (coss*NV_Ith_S(y, IcIm) - sinn*NV_Ith_S(y, IcRe));	// b1 Re
     NV_Ith_S(ydot, IbIm) -= (sinn*NV_Ith_S(y, IcIm) + coss*NV_Ith_S(y, IcRe));	// b1 Im
#ifdef DEBUGf
     cout << endl << "IcRe " << IcRe << " IcIm " << IcIm << " IbRe " << IbRe << " IbIm ";
     cout << IbIm << " Vcb " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
#endif
    }
  for (i = 0; i < Nb - 1; i++)				// Vb_nb_{n+1}
   for (n = 0; n < N_vib; n++)
    for (m = 0; m < N_vib; m++) {
     IbRe = Ib_vib + i*N_vib + n;
     IbIm = IbRe + NEQ_vib;
     IBRe = Ib_vib + (i+1)*N_vib + m;
     IBIm = IBRe + NEQ_vib;
     Vee = Vbridge[i+1]*FCbb[n][m];
     coss = Vee*cos((energy[IbRe] - energy[IBRe])*t);
     sinn = Vee*sin((energy[IbRe] - energy[IBRe])*t);
     NV_Ith_S(ydot, IbRe) += (coss*NV_Ith_S(y, IBIm) + sinn*NV_Ith_S(y, IBRe));	// b Re
     NV_Ith_S(ydot, IbIm) += (sinn*NV_Ith_S(y, IBIm) - coss*NV_Ith_S(y, IBRe));	// b Im
     NV_Ith_S(ydot, IBRe) += (coss*NV_Ith_S(y, IbIm) - sinn*NV_Ith_S(y, IbRe));	// b1 Re
     NV_Ith_S(ydot, IBIm) -= (sinn*NV_Ith_S(y, IbIm) + coss*NV_Ith_S(y, IbRe));	// b1 Im
#ifdef DEBUGf
     cout << endl << "IbRe " << IbRe << " IbIm " << IbIm << " IBRe " << IBRe << " IBIm ";
     cout << IBIm << " Vbb " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
#endif
    }
 }
 else {						// no bridge
  realtype Vkc = V[Ik][Ic];
  for (i = 0; i < Nk; i++)
   for (j = 0; j < Nc; j++)
    for (n = 0; n < N_vib; n++)
     for (m = 0; m < N_vib; m++) {
      IkRe = Ik_vib + i*N_vib + n;			// indices
      IkIm = IkRe + NEQ_vib;
      IcRe = Ic_vib + j*N_vib + m;
      IcIm = IcRe + NEQ_vib;
      /* Franck-Condon indices should be consistent in direction from one end of the 
       * system to the other.  That is, if the first index is for the beginning state
       * the second should be for the next in the chain, and so on. */
      Vee = Vkc*FCkc[n][m];
      coss = Vee*cos((energy[IkRe] - energy[IcRe])*t);	// do the cosine now for speed
      sinn = Vee*sin((energy[IkRe] - energy[IcRe])*t);	// do the sine now for speed
      NV_Ith_S(ydot, IkRe) += (coss*NV_Ith_S(y, IcIm) + sinn*NV_Ith_S(y, IcRe)); // k Re
      NV_Ith_S(ydot, IkIm) += (sinn*NV_Ith_S(y, IcIm) - coss*NV_Ith_S(y, IcRe)); // k Im
      NV_Ith_S(ydot, IcRe) += (coss*NV_Ith_S(y, IkIm) - sinn*NV_Ith_S(y, IkRe)); // c Re
      NV_Ith_S(ydot, IcIm) -= (sinn*NV_Ith_S(y, IkIm) + coss*NV_Ith_S(y, IkRe)); // c Im
#ifdef DEBUGf
      cout << endl << "IkRe " << IkRe << " IkIm " << IkIm << " IcRe " << IcRe << " IcIm ";
      cout << IcIm << " V " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
#endif
     }
 }

#ifdef DEBUGf
 cout << "\n\nN_Vectors at time " << t << ":\n";
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Re[k(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ik_vib+i*N_vib+j);
   cout << ", ydot = " << NV_Ith_S(ydot, Ik_vib+i*N_vib+j) << endl;
  }
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Re[c(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ic_vib+i*N_vib+j);
   cout << ", ydot = " << NV_Ith_S(ydot, Ic_vib+i*N_vib+j) << endl;
  }
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Re[b(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ib_vib+i*N_vib+j);
   cout << ", ydot = " << NV_Ith_S(ydot, Ib_vib+i*N_vib+j) << endl;
  }
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Im[k(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ik_vib+i*N_vib+j+NEQ_vib);
   cout << ", ydot = " << NV_Ith_S(ydot, Ik_vib+i*N_vib+j+NEQ_vib) << endl;
  }
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Im[c(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ic_vib+i*N_vib+j+NEQ_vib);
   cout << ", ydot = " << NV_Ith_S(ydot, Ic_vib+i*N_vib+j+NEQ_vib) << endl;
  }
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++) {
   cout << "Im[b(" << i << "," << j << ")]: y = " << NV_Ith_S(y, Ib_vib+i*N_vib+j+NEQ_vib);
   cout << ", ydot = " << NV_Ith_S(ydot, Ib_vib+i*N_vib+j+NEQ_vib) << endl;
  }
#endif

#ifdef DEBUG_SAI
   // this bit spits out dy/dt for t = 0, so one can check initial conditions.
   // the file dydt.out can be compared directly to dydt.dat from Sai's code
 double end_time = 3.445;
 if (t == 0 || (t > 0.0001 && t < end_time && t != last_t)) {
  cout << "tee is " << t << endl;
  FILE * dydt;
  dydt = fopen("dydt.out", "a+");
  cout << "whoot\n";
  // for (i = 0; i < 2*NEQ_vib; i++)
   // cout << NV_Ith_S(ydot, i)*41.3414 << endl;
  for (i = 0; i < Nc ; i++)
   for (j = 0; j < N_vib ; j++) {
    fprintf(dydt, "Re(c%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ic_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(c%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ic_vib + i*N_vib + j + NEQ_vib)*41.3414);
   }
  for (i = 0; i < Nb ; i++)
   for (j = 0; j < N_vib ; j++) {
    fprintf(dydt, "Re(b%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ib_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(b%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ib_vib + i*N_vib + j + NEQ_vib)*41.3414);
   }
  for (i = 0; i < Nk ; i++)
   for (j = 0; j < N_vib ; j++) {
    fprintf(dydt, "Re(k%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ik_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(k%d,%-.5g): %-.9g\n", j, t, NV_Ith_S(ydot, Ik_vib + i*N_vib + j + NEQ_vib)*41.3414);
   }
  fclose(dydt);
  last_t = t;
 }
 // this is when I'm looking for a sign error somewhere... *sign*
 if (t == -1) {
  NV_Ith_S(ydot, 7) = -NV_Ith_S(ydot, 7);
  cout << NV_Ith_S(ydot, 7)*41.3414 << endl;
 }
#endif

 return 0;
}


int Normalize_NV(N_Vector nv, realtype total) {
// normalizes an N_Vector so that the populations of all states are
// normalized to value 'total'

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


int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep) {
// compute the six-point derivative of an array.
// Assumes array elements are evenly spaced.
// Output array is six elements shorter than input.

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


int Output_checkpoint(
#ifdef DEBUG
  FILE * realImaginary, 
#endif
  double ** allprobs, N_Vector outputData, realtype time,
  realtype * totK, realtype * totL, realtype * totC, realtype * totB, realtype ** vibProb, realtype * times,
  realtype * qd_est, realtype * qd_est_diag, realtype * energy_expectation, int index, realtype * energies,
  double kBandEdge, double kBandTop, double * k_pops) {
// computes many time-dependent properties

 int i, j, k, l;				// counters
 int Re1, Im1, Re2, Im2;			// indices for real and imaginary components
 double sinn, coss;				// these are used for computing observables
 						// in the interaction picture.
 int Idx;
 double sumkpop = 0;
 double sumcpop = 0;
 double sumbpop = 0;
 double sumlpop = 0;
 double temp;

#ifdef DEBUG
 for (i = 0; i < NEQ_vib; i++) {
  fprintf(realImaginary, "%-.9e %-.9e\n", NV_Ith_S(outputData, i), NV_Ith_S(outputData, i+NEQ_vib));
 }
 fprintf(realImaginary, "\n");
#endif

 for (i = 0; i < Nk; i++) {			// k populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ik_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  allprobs[index][Ik+i] = temp;
  sumkpop += temp;
 }
 for (i = 0; i < Nl; i++) {			// l populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Il_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  allprobs[index][Il+i] = temp;
  sumlpop += temp;
 }
 for (i = 0; i < Nc; i++) {			// c populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ic_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  allprobs[index][Ic+i] = temp;
  sumcpop += temp;
 }
 for (i = 0; i < Nb; i++) {			// b populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ib_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  sumbpop += temp;
  allprobs[index][Ib+i] = temp;
 }
 for (i = 0; i < N_vib; i++) {
  temp = 0;
  for (j = 0; j < NEQ; j++) {
   temp += pow(NV_Ith_S(outputData,i+j*N_vib),2) + pow(NV_Ith_S(outputData, i+j*N_vib+NEQ_vib),2);
  }
  vibProb[index][i] = temp;
 }
 totK[index] = sumkpop;
 totL[index] = sumlpop;
 totC[index] = sumcpop;
 totB[index] = sumbpop;
 times[index] = time;

 temp = 0.0;
 for (i = 0; i < NEQ; i++) {		// loop over all states
  for (j = 0; j < N_vib; j++) {		// loop over all vibronic states
   Re1 = i*N_vib + j;			// index for current state
   Im1 = i*N_vib + j + NEQ_vib;	
   temp += energy[Re1]*(pow(NV_Ith_S(outputData,Re1),2) + pow(NV_Ith_S(outputData,Im1),2));
   for (k = 0; k < NEQ; k++) {		// loop over all states (for coupling)
    for (l = 0; l < N_vib; l++) {	// loop over all vibronic states (for coupling)
     Re2 = k*N_vib + l;			// index for coupled state
     Im2 = k*N_vib + l + NEQ_vib;
     // WARNING: this only works correctly when the Franck-Condon displacement
     // factors are all 0. In order for this to be corrected, I will have to
     // rewrite this section to account for the differing FC factors for
     // different states (k, c, b...)
     sinn = V[i][k]*sin((energy[Re1] - energy[Re2])*time);	// precalculate sin and cos to save space
     coss = V[i][k]*cos((energy[Re1] - energy[Re2])*time);
     temp += NV_Ith_S(outputData,Re1)*NV_Ith_S(outputData,Re2)*coss;
     temp += NV_Ith_S(outputData,Im1)*NV_Ith_S(outputData,Im2)*coss;
     temp -= NV_Ith_S(outputData,Re1)*NV_Ith_S(outputData,Im2)*sinn;
     temp += NV_Ith_S(outputData,Im1)*NV_Ith_S(outputData,Re2)*sinn;
    }
   }
  }
 }
 energy_expectation[index] = temp;

 return 0;
}


int Analytical_c (
 double tout, int timesteps, realtype * energies,
 double kBandEdge, double kBandTop, N_Vector y) {
// computes analytically the population on a single c state

 // energy spacing in bulk
 complex <double> dE ((kBandTop-kBandEdge)/(Nk-1), 0);
 // bulk-QD coupling
 complex <double> Vee (V[Ik][Ic], 0);
 // rate constant (can be defined also as K/2)
 complex <double> K = complex <double> (3.1415926535,0)*pow(Vee,2)/dE;
 // time
 complex <double> t (0, 0);
 // energy differences
 complex <double> wnm (0, 0);
 complex <double> wnnp (0, 0);
 complex <double> wnpm (0, 0);
 // coefficients
 complex <double> cm (0, 0);
 complex <double> cn (0, 0);
 complex <double> cn_term1 (0, 0);
 complex <double> cn_term2 (0, 0);
 complex <double> cn_diag (0, 0);
 complex <double> cn_offdiag (0, 0);
 double cn_tot;
 // complex numbers are dumb
 complex <double> II (0, 1);
 complex <double> I2 (-1, 0);

 FILE * ms_est;
 FILE * ms_est_tot;
 FILE * c_diag;			// diagonal portion of coefficients
 FILE * c_offdiag;		// diagonal portion of coefficients

 ms_est = fopen("ms_est.out", "w");
 ms_est_tot = fopen("ms_est_tot.out", "w");
 c_diag = fopen("c_diag.out", "w");
 c_offdiag = fopen("c_offdiag.out", "w");

 for (t = complex<double>(0,0); real(t) <= tout; t += complex<double>(tout/timesteps,0)) {
  cn_tot = 0.0;
  fprintf(ms_est, "%.7g", real(t));
  fprintf(c_diag, "%.7g", real(t));
  fprintf(c_offdiag, "%.7g", real(t));
  for (int n = 0; n < Nc; n++) {
   cn = complex<double>(0, 0);
   cn_diag = complex<double>(0, 0);
   cn_offdiag = complex<double>(0, 0);
   cn_term1 = complex<double>(0, 0);
   cn_term2 = complex<double>(0, 0);
   for (int m = 0; m < Nk; m++) {
    cm = complex<double>(NV_Ith_S(y,m), 0);
    wnm = complex<double>(energies[Ic + n] - energies[Ik + m], 0);
    cn -= II*Vee*sqrt(dE)*cm // first term
       *(exp(II*wnm*t) - exp(I2*K*t))
       /(K + II*wnm);
    // eh, this may not be the way to look at the (off-)diagonal terms
    cn_term1 -= II*Vee*sqrt(dE)*cm*(exp(II*wnm*t))/(K + II*wnm); // 1st part of first term
    cn_term2 += II*Vee*sqrt(dE)*cm*(exp(I2*K*t))/(K + II*wnm); // 2nd part of first term
    cn_diag += cn_term1*conj(cn_term1) + cn_term2*conj(cn_term2);
    cn_offdiag += cn_term1*conj(cn_term2) + cn_term2*conj(cn_term1);
    for (int np = 0; np < Nc; np++) {
     if (np == n) continue;
     wnnp = complex<double>(energies[Ic + n] - energies[Ic + np], 0);
     wnpm = complex<double>(energies[Ic + np] - energies[Ik + m], 0);
     cn += II*Vee*K*sqrt(dE)*cm
	*((exp(II*wnm*t) - exp(I2*K*t))
	  /((K+II*wnpm)*(K + II*wnm))
	  -(exp(I2*(K-II*wnnp)*t) - exp(I2*K*t))
	  /((K+II*wnpm)*II*wnnp));
    }
   }

   complex<double>thexfactor (1.0/pow(sqrt((kBandTop-kBandEdge)/(Nk-1)),2), 0);
   // complex<double>thexfactor (1, 0);
   fprintf(ms_est, " %.7g", real(thexfactor*cn*conj(cn)));
   fprintf(c_diag, " %.7g", real(thexfactor*cn_diag));
   fprintf(c_offdiag, " %.7g", real(thexfactor*cn_offdiag));
   cn_tot += real(thexfactor*cn*conj(cn));
  }
  fprintf(ms_est, "\n");
  fprintf(c_diag, "\n");
  fprintf(c_offdiag, "\n");
  fprintf(ms_est_tot, "%.7g %.7g\n", real(t), cn_tot);
 }

 fclose(ms_est);
 fclose(ms_est_tot);
 fclose(c_diag);
 fclose(c_offdiag);

 return 0;
}
realtype Integrate_arrays (realtype * values, realtype * time, int num) {
// Riemann sum of an array (values) at time points (time).
// Does not assume equal spacing in time.

 int i;
 realtype riemann = 0;

 for (i = 0; i < num-1; i++)
  riemann += (values[i+1] + values[i])*(time[i+1]-time[i])/2;

 return riemann;
}


realtype Find_array_maximum (realtype * inputArray, int num) {
// Returns maximum element in an array.

 int i;
 realtype currentMax = inputArray[0];

 for (i = 1; i < num; i++)
  if (inputArray[i] > currentMax)
   currentMax = inputArray[i];

 return currentMax;
}


realtype Find_first_array_maximum (realtype * inputArray, int num) {
// Finds the first maximum in an array (the first point where the next
// point is smaller in value).

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


int Find_first_array_maximum_index (realtype * inputArray, int num) {
// This function returns the index of the first maximum in an array.
// Warning: will return 0 as index of first maximum if the second array element
// is less than the first.  This may not be what you want.

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


int Find_array_maximum_index (realtype * inputArray, int num) {
// returns index of first maximum in an array.

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


void Compute_final_outputs (double ** allprobs, realtype * time, realtype * tk,
  realtype * tl, realtype * tc, realtype * tb, realtype ** vibProb, realtype * energies,
  realtype * energy_expectation, int num, double * qd_est, double * qd_est_diag) {
// makes output files for time-dependent properties

 FILE * tkprob;
 FILE * tlprob;
 FILE * tcprob;
 FILE * vibprob;
 FILE * kprobs;
 FILE * kprobs_gnuplot;
 FILE * cprobs;
 FILE * cprobs_gnuplot;
 FILE * bprobs;
 FILE * Ikprob;
 FILE * Icprob;
 FILE * kmax;
 FILE * cmax;
 FILE * cmax_t;
 FILE * cmax_first;
 FILE * cmax_first_t;
 FILE * totprob;
 FILE * energy;
 FILE * times;
 FILE * pump_intensity;
 FILE * energy_exp;
 FILE * qd_estimate;
 FILE * qd_estimate_diag;
 int i, j;
 realtype summ;

 tkprob = fopen("tkprob.out", "w");
 tlprob = fopen("tlprob.out", "w");
 tcprob = fopen("tcprob.out", "w");
 vibprob = fopen("vibprob.out", "w");
 kprobs = fopen("kprobs.out", "w");
 kprobs_gnuplot = fopen("kprobs_gnuplot.out", "w");
 cprobs_gnuplot = fopen("cprobs_gnuplot.out", "w");
 cprobs = fopen("cprobs.out", "w");
 bprobs = fopen("bprobs.out", "w");
 totprob = fopen("totprob.out", "w");
 Ikprob = fopen("Ikprob.out", "w");
 Icprob = fopen("Icprob.out", "w");
 kmax = fopen("kmax.out", "w");
 cmax = fopen("cmax.out", "w");
 cmax_t = fopen("cmax_t.out", "w");
 cmax_first = fopen("cmax_first.out", "w");
 cmax_first_t = fopen("cmax_first_t.out", "w");
 energy = fopen("energy.out", "w");
 times = fopen("times.out", "w");
 pump_intensity = fopen("pump_intensity.out", "w");
 energy_exp = fopen("energy_exp.out", "w");
 qd_estimate = fopen("qd_est.out", "w");
 qd_estimate_diag = fopen("qd_est_diag.out", "w");

 for (i = 0 ; i < num ; i++) {				// print k probabilities over time
  fprintf(kprobs, "%-.7g", time[i]);
  for (j = 0; j < Nk ; j++)
   fprintf(kprobs, " %-.7g", allprobs[i][Ik+j]);
  fprintf(kprobs, "\n");
 }
 for (i = 0 ; i < num ; i++) {
  for (j = 0 ; j < Nk ; j++ )
   fprintf(kprobs_gnuplot, "%-.7g %-.7g %-.7g\n", time[i], energies[Ik_vib + j*N_vib], allprobs[i][Ik+j]);
  if (i < (num - 1))
   fprintf(kprobs_gnuplot, "\n");			// makes a blank line for gnuplot
 }


 for (i = 0 ; i < num ; i++) {				// print c probabilities over time
  fprintf(cprobs, "%-.7g", time[i]);
  for (j = 0; j < Nc ; j++)
   fprintf(cprobs, " %-.7g", allprobs[i][Ic+j]);
  fprintf(cprobs, "\n");
 }
 for (i = 0 ; i < num ; i++) {
  for (j = 0 ; j < Nc ; j++ )
   fprintf(cprobs_gnuplot, "%-.7g %-.7g %-.7g\n", time[i], energies[Ic_vib + j*N_vib], allprobs[i][Ic+j]);
  if (i < (num - 1))
   fprintf(cprobs_gnuplot, "\n");			// makes a blank line for gnuplot
 }

 for (i = 0 ; i < num ; i++) {				// print b probabilities over time
  fprintf(bprobs, "%-.7g", time[i]);
  for (j = 0; j < Nb ; j++)
   fprintf(bprobs, " %-.7g", allprobs[i][Ib+j]);
  fprintf(bprobs, "\n");
 }

 if (Nb > 0) {
  FILE * tbprob;
  FILE * Ibprob;
  FILE * bmax;
  double max_b_prob = 0;

  tbprob = fopen("tbprob.out", "w");
  Ibprob = fopen("Ibprob.out", "w");
  bmax = fopen("bmax.out", "w");

  for (i = 0; i <= num; i++)				// print total b population
   fprintf(tbprob, "%-.7g %-.7g\n", time[i], tb[i]);
  
  fprintf(Ibprob, "%-.7g", Integrate_arrays(tb, time, num+1));
  
  for (i = 0 ; i < num + 1 ; i++) {
   for (j = 0 ; j < Nb ; j++) {
    if (allprobs[i][Ib+j] > max_b_prob) {
     max_b_prob = allprobs[i][Ib+j];
    }
   }
  }
  //fprintf(bmax, "%-.7g", Find_array_maximum(tb, num+1));
  fprintf(bmax, "%-.7g", max_b_prob);
  
  fclose(tbprob);
  fclose(Ibprob);
  fclose(bmax);
 }

 for (i = 0; i <= num; i++) {				// print total k, c population
  fprintf(tkprob, "%-.7g %-.7g\n", time[i], tk[i]);
  fprintf(tlprob, "%-.7g %-.7g\n", time[i], tl[i]);
  fprintf(tcprob, "%-.7g %-.7g\n", time[i], tc[i]);
  summ = tk[i] + tc[i] + tl[i];
  if (Nb > 0)
   summ += tb[i];
  fprintf(totprob, "%-.7g %-.15g\n", time[i], summ);
 }

 for (i = 0; i <= num; i++) {				// print vibrational populations
  fprintf(vibprob, "%-.7g %-.7g", time[i], vibProb[i][0]);
  for (j = 1; j < N_vib; j++)
   fprintf(vibprob, " %-.7g", vibProb[i][j]);
  fprintf(vibprob, "\n");
 }

 fprintf(Ikprob, "%-.7g", Integrate_arrays(tk, time, num+1));
 fprintf(Icprob, "%-.7g", Integrate_arrays(tc, time, num+1));
 
 fprintf(kmax, "%-.7g", Find_array_maximum(tk, num+1));
 fprintf(cmax, "%-.7g", Find_array_maximum(tc, num+1));
 fprintf(cmax_first, "%-.7g", Find_first_array_maximum(tc, num+1));
 fprintf(cmax_t, "%-.7g", time[Find_array_maximum_index(tc, num+1)]);
 fprintf(cmax_first_t, "%-.7g", time[Find_first_array_maximum_index(tc, num+1)]);

 // energy.out should be all the energies on one row, since it's used for
 // the movie maker.
 fprintf(energy, "%-.7g", energies[0]);
 for (i = 1; i < NEQ; i++)
  fprintf(energy, " %-.7g", energies[i*N_vib]);

 for (i = 0; i <= num; i++) {
  fprintf(times, "%-.7g\n", time[i]);
  fprintf(pump_intensity, "%-.7g %-.7g\n", time[i], pump(time[i]));
  fprintf(energy_exp, "%-.7g %-.7g\n", time[i], energy_expectation[i]);
  fprintf(qd_estimate, "%-.7g %-.7g\n", time[i], qd_est[i]);
  fprintf(qd_estimate_diag, "%-.7g %-.7g\n", time[i], qd_est_diag[i]);
 }

 fclose(tkprob);
 fclose(tlprob);
 fclose(tcprob);
 fclose(vibprob);
 fclose(kprobs);
 fclose(kprobs_gnuplot);
 fclose(cprobs);
 fclose(cprobs_gnuplot);
 fclose(bprobs);
 fclose(totprob);
 fclose(Ikprob);
 fclose(Icprob);
 fclose(kmax);
 fclose(cmax);
 fclose(cmax_first);
 fclose(cmax_t);
 fclose(cmax_first_t);
 fclose(energy);
 fclose(times);
 fclose(qd_estimate);
 fclose(qd_estimate_diag);

 // compute derivatives
 FILE * tkDeriv;
 FILE * tcDeriv;
 double * tkderivs = new double [numOutputSteps-5];
 double * tcderivs = new double [numOutputSteps-5];
 tkDeriv = fopen("tkderiv.out","w");
 tcDeriv = fopen("tcderiv.out","w");
 Derivative(tk, numOutputSteps, tkderivs, time[1]-time[0]);
 Derivative(tc, numOutputSteps, tcderivs, time[1]-time[0]);
 for (i = 0; i < numOutputSteps-5; i++) {
  fprintf(tkDeriv, "%-.7g %-.9g\n", time[i+2], tkderivs[i]);
  fprintf(tcDeriv, "%-.7g %-.9g\n", time[i+2], tcderivs[i]);
 }
 fclose(tkDeriv);
 fclose(tcDeriv);
 delete [] tkderivs;
 delete [] tcderivs;

 // compute rates.  Rate is calculated as (dP(t)/dt)/P(t), where P(t) is
 // the population at time t.
 FILE * tkRate;
 FILE * tcRate;
 tkRate = fopen("tkrate.out","w");
 tcRate = fopen("tcrate.out","w");
 for (i = 0; i < numOutputSteps-5; i++) {
  fprintf(tkRate, "%-.7g %-.9g\n", time[i+2], tkderivs[i]/tk[i+2]);
  fprintf(tcRate, "%-.7g %-.9g\n", time[i+2], tcderivs[i]/tc[i+2]);
 }
 fclose(tkRate);
 fclose(tcRate);

}


void outputYData(realtype * ydata, int n) {
// prints out initial wave function.  Inputs are the wave function array and
// the number of equations.
 FILE * psi_start;
 psi_start = fopen("psi_start.out", "w");
 for (int i = 0; i < n; i++) {
  fprintf(psi_start, "%-.9g %-.9g\n", ydata[i],2, ydata[i+n]);
 }
 fclose(psi_start);
}

void buildHamiltonian(realtype * H, realtype * energy, realtype ** V, int N, int N_vib,
                      realtype ** FCkb, realtype ** FCbb, realtype ** FCbc) {
// builds a Hamiltonian from site energies and couplings.
 int i, j;	// counters!
 int n, m;	// more counters!
 int idx1, idx2;
 
 if (bridge_on) {
  // assign diagonal elements
#ifdef DEBUG
  fprintf(stderr, "Assigning diagonal elements in Hamiltonian (with bridge).\n");
#endif
  for (i = 0; i < NEQ_vib; i++) {
   H[i*N + i] = energy[i];
  }
  // assign bulk-bridge coupling
  // General way to find index in vibrational matrix:
  // H_{a_i,v_n;b_j,u_m} (ith state of type a, vib n; jth state of type b, vib m)
  // (idx1 + n)*N + idx2 + m
  // ((Ia_vib + i)*N_vib + n)*N + (Ib_vib + j)*N_vib + m
#ifdef DEBUG
  fprintf(stderr, "Assigning bulk-bridge coupling elements in Hamiltonian.\n");
#endif
  idx2 = Ib_vib;
  for (i = 0; i < Nk; i++) {
   idx1 = (Ik + i)*N_vib;
   for (n = 0; n < N_vib; n++) {
    for (m = 0; m < N_vib; m++) {
     //H[((Ik_vib + i)*N_vib + n)*N + (Ib_vib + m)] = V[(Ik_vib + i)*N_vib + n][Ib_vib + m]*FCkb[n][m];
#ifdef DEBUG2
     fprintf(stderr, "H[(%d + %d)*%d + %d + %d] = ", idx1, n, N, idx2, m);
     fprintf(stderr, "V[%d][%d]*FCkb[%d][%d] = ", idx1, idx2, n, m);
     fprintf(stderr, "%e\n", V[Ik+i][Ib]*FCkb[n][m]);
#endif
     H[(idx1 + n)*N + idx2 + m] = V[Ik+i][Ib]*FCkb[m][n];
     H[(idx2 + m)*N + idx1 + n] = V[Ib][Ik+i]*FCkb[n][m];
    }
   }
  }
  // assign bridge-bridge couplings
#ifdef DEBUG
  fprintf(stderr, "Assigning bridge-bridge coupling elements in Hamiltonian.\n");
#endif
  for (i = 1; i < Nb; i++) {
   idx1 = (Ib + i)*N_vib;
   idx2 = (Ib+ i + 1)*N_vib;
   for (n = 0; n < N_vib; n++) {
    for (m = 0; m < N_vib; m++) {
#ifdef DEBUG2
     fprintf(stderr, "H[(%d + %d)*%d + %d + %d] = ", idx1, n, N, idx2, m);
     fprintf(stderr, "V[%d][%d]*FCkb[%d][%d] = ", idx1, idx2, n, m);
     fprintf(stderr, "%e\n", V[Ib+i][Ib+i+1]*FCbb[n][m]);
#endif
     H[(idx1 + n)*N + idx2 + m] = V[Ib+i][Ib+i+1]*FCbb[n][m];
     H[(idx2 + m)*N + idx1 + n] = V[Ib+i+1][Ib+i]*FCbb[m][n];
    }
   }
  }
  // assign bridge-QD coupling
#ifdef DEBUG
  fprintf(stderr, "Assigning bridge-QD coupling elements in Hamiltonian.\n");
#endif
  idx2 = (Ib + Nb - 1)*N_vib;
  for (i = 0; i < Nc; i++) {
   idx1 = (Ic + i)*N_vib;
   for (n = 0; n < N_vib; n++) {
    for (m = 0; m < N_vib; m++) {
#ifdef DEBUG2
     fprintf(stderr, "H[(%d + %d)*%d + %d + %d] = ", idx1, n, N, idx2, m);
     fprintf(stderr, "V[%d][%d]*FCbc[%d][%d] = ", Ic+i, Ib+Nb-1, n, m);
     fprintf(stderr, "%e\n", V[Ic+i][Ib+Nb-1]*FCbc[n][m]);
#endif
     H[(idx1 + n)*N + idx2 + m] = V[Ic+i][Ib+Nb-1]*FCbc[n][m];
     H[(idx2 + m)*N + idx1 + n] = V[Ib+Nb-1][Ic+i]*FCbc[m][n];
    }
   }
  }
 }
 else {
#ifdef DEBUG
  fprintf(stderr, "Assigning elements in Hamiltonian (no bridge).\n");
#endif
  for (i = 0; i < N; i++) {
   // diagonal
   H[i*N + i] = energy[i];
#ifdef DEBUG
   cout << "diagonal element " << i << " of H is " << energy[i] << "\n";
#endif
   for (j = 0; j < i; j++) {
    // off-diagonal
    H[i*N + j] = V[i][j];
    H[j*N + i] = V[j][i];
   }
  }
 }
}

void printVector(realtype * W, int N, char * fileName) {
// prints a vector W of length N
 int i;        // counter
 FILE * out;   // output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e\n", W[i]);
 }

 fclose(out);
}

void printPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName) {
// prints the elementwise square of a complex vector
 int i;		// counter!
 FILE * out;	// output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e %-.9e\n", evals[i], (pow(v[i].re,2) + pow(v[i].im,2)));
 }
 
 fclose(out); 
}

void printCVector(complex16 * v, int N, char * fileName) {
// prints a complex vector v of length N
 int i;		// counter!
 FILE * out;	// output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e %-.9e\n", v[i].re, v[i].im);
 }

 fclose(out);
}

void printCVectorTime(complex16 * v, int N, int M, char * fileName) {
// prints a complex vector v of length N with M time steps
 int i, j;	// counters!
 FILE * out;	// output file

 out = fopen(fileName, "w");

 for (j = 0; j < M; j++) {
  for (i = 0; i < N; i++) {
   fprintf(out, "%-.9e %-.9e\n", v[j*N+i].re, v[j*N+i].im);
  }
  fprintf(out, "\n");
 }

 fclose(out);
}

void printSquareMatrix(realtype * M, int N, char * fileName) {
// prints a square matrix M of dimension N
 int i, j;     // counters
 FILE * out; // output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  for (j = 0; j < N; j++) {
   if (j == 0) {
    fprintf(out, "%-.9e", M[i*N + j]);
   }
   else {
    fprintf(out, " %-.9e", M[i*N + j]);
   }
  }
  fprintf(out, "\n");
 }

 fclose(out);
}

void projectSubsystems(realtype * evecs, realtype * evals, int dim) {
// projects the initial wavefunction from the state basis onto the
// subsystems (bulk, bridge, QD, etc.)
 int i, j;	// counters!
 // create arrays for bulk, bridge, QD
 double * bu_proj = new double [dim];
 double * br_proj = new double [dim];
 double * qd_proj = new double [dim];
 // sum up subsystem projections for each eigenvector
 for (i = 0; i < dim; i++) {
  bu_proj[i] = 0.0;
  for (j = Ik; j < Ik + Nk; j++) {
   bu_proj[i] += pow(evecs[i*dim + j],2);
  }
  br_proj[i] = 0.0;
  for (j = Ib; j < Ib + Nb; j++) {
   br_proj[i] += pow(evecs[i*dim + j],2);
  }
  qd_proj[i] = 0.0;
  for (j = Ic; j < Ic + Nc; j++) {
   qd_proj[i] += pow(evecs[i*dim + j],2);
  }
 }

 // write projections to file
 FILE * projections;
 projections = fopen("projections.out", "w");
 for (i = 0; i < dim; i++) {
  fprintf(projections, "%-.9g %-.9g %-.9g %-.9g\n",
          evals[i], bu_proj[i], br_proj[i], qd_proj[i]);
 }
 fclose(projections);

 // make a plotting file
 FILE * proj_plot;
 proj_plot = fopen("projections.plt", "w");
 fprintf(proj_plot, "#!/usr/bin/env gnuplot\n");
 fprintf(proj_plot, "set terminal postscript enhanced color size 20cm,14cm font 'Courier-Bold,12'\n");
 fprintf(proj_plot, "set output 'pdos_stack.eps'\n");
 fprintf(proj_plot, "set key outside right\n");
 fprintf(proj_plot, "set tics scale 0\n");
 fprintf(proj_plot, "set yr [0:1]\n");
 fprintf(proj_plot, "set ytics ('0' 0, '1' 1)\n");
 fprintf(proj_plot, "set rmargin 18\n");
 fprintf(proj_plot, "set multiplot layout 4,1\n");
 fprintf(proj_plot, "plot '../outs/projections.out' u 1:2 w im lw 2 lc rgb 'red' t 'Bulk'\n");
 fprintf(proj_plot, "plot '../outs/projections.out' u 1:3 w im lw 2 lc rgb 'blue' t 'Bridge'\n");
 fprintf(proj_plot, "plot '../outs/projections.out' u 1:4 w im lw 2 lc rgb 'green' t 'QD'\n");
 fprintf(proj_plot, "set xlabel 'Energy (E_h)'\n");
 fprintf(proj_plot, "plot '../outs/psi2_start_e.out' w im lw 2 lc rgb 'black' t 'Psi'\n");
 fprintf(proj_plot, "unset multiplot\n");
 fclose(proj_plot);

 // clean up
 delete [] bu_proj;
 delete [] br_proj;
 delete [] qd_proj;
}

void projectStateToSite(complex16 * psi_E_t, int dim, realtype * evecs, complex16 * psi_S_t, int timesteps) {
// projects the time-dependent wavefunction from the state basis to the site basis
 int i;
 int M = dim;
 int N = timesteps+1;
 int K = dim;
 complex16 ALPHA;
 ALPHA.re = 1.0;
 ALPHA.im = 0.0;
 // convert evecs to a complex matrix A
 complex16 * A = new complex16 [dim*dim];
 for (i = 0; i < dim*dim; i++) {
  A[i].re = evecs[i];
  A[i].im = 0.0;
 }
 int LDA = dim;
 int LDB = dim;
 complex16 BETA;
 BETA.re = 0.0;
 BETA.im = 0.0;
 int LDC = dim;
 int INCX = 1;
 int INCY = 1;
 // designating evecs col-major, essentially calling it the transpose
 cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, K, &ALPHA,
  A, LDA, psi_E_t, LDB, &BETA, psi_S_t, LDC);
 // can also go row-by-row.
 //for (i = 0; i <=timesteps; i++) {
 // cblas_zgemv(CblasRowMajor, CblasTrans, M, N, &ALPHA, A, LDA, psi_E_t+(i*dim),
 //             INCX, &BETA, psi_S_t+(i*dim), INCY);
 //}
 delete [] A;
}

void projectSiteToState(complex16 * psi_S, int dim, realtype * evecs, complex16 * psi_E) {
// projects the wavefunction in the site basis (psi_S) onto the eigenstate
// basis (psi_E).  Requires the orthonormal eigenvectors (evecs), and assumes
// they are in a column matrix format (row-major)
 int i;			// counter!
 // BLAS variables
 int M = dim;
 int N = dim;
 complex16 ALPHA;
 ALPHA.re = 1.0;
 ALPHA.im = 0.0;
 // convert the real matrix of evecs to a complex matrix A
 complex16 * A = new complex16 [dim*dim];
 for (i = 0; i < dim*dim; i++) {
  A[i].re = evecs[i];
  A[i].im = 0.0;
 }
 int LDA = dim;
 int INCX = 1;
 complex16 BETA;
 BETA.re = 0.0;
 BETA.im = 0.0;
 int INCY = 1;

 cblas_zgemv(CblasRowMajor, CblasNoTrans, M, N, &ALPHA, A, LDA, psi_S, INCX,
  &BETA, psi_E, INCY);
 delete [] A;
}

void propagatePsi(complex16 * psi_E, complex16 * psi_E_t, int N,
 double * evals, int timesteps, double tout) {
// propagates a wavefunction in time in the eigenbasis (really just applying
// a phase factor).
 int i, j;	// counters!
 double t;	// time for each step

 // get t = 0 for free
 for (i = 0; i < N; i++) {
  psi_E_t[i].re = psi_E[i].re;
  psi_E_t[i].im = psi_E[i].im;
 }
 // get the rest of the timesteps
 for (j = 1; j <= timesteps; j++) {
  t = tout*(((double)j)/((double)timesteps));
  for (i = 0; i < N; i++) {
#ifdef DEBUGf
   if (i == 0) {
    cout << "timestep number " << j << " is " << t << "\n";
   }
   cout << "assigning psi_E_t[" << j << "," << i << "]...\n";
#endif
   psi_E_t[j*N + i].re = psi_E[i].re*cos(evals[i]*t) - psi_E[i].im*sin(evals[i]*t);
   psi_E_t[j*N + i].im = psi_E[i].im*cos(evals[i]*t) + psi_E[i].re*sin(evals[i]*t);
  }
 }
}

void makeOutputsTI(complex16 * psi_t, int dim, double * t, int timesteps) {
// makes outputs for the result of a time-independent H time propagation
 FILE * tkprob;
 FILE * tcprob;
 FILE * tbprob;
 FILE * tlprob;
 FILE * kprobs;
 FILE * cprobs;
 FILE * bprobs;
 FILE * lprobs;
 FILE * kmax;
 FILE * cmax;
 FILE * cmax_t;
 FILE * cmax_first;
 FILE * cmax_first_t;
 FILE * totprob;
 FILE * energy;
 FILE * times;
 int i, j;
 double summ;
 double maxx;
 int maxx_t;

 tkprob = fopen("tkprob.out", "w");
 tcprob = fopen("tcprob.out", "w");
 tbprob = fopen("tbprob.out", "w");
 tlprob = fopen("tlprob.out", "w");
 kprobs = fopen("kprobs.out", "w");
 cprobs = fopen("cprobs.out", "w");
 bprobs = fopen("bprobs.out", "w");
 lprobs = fopen("lprobs.out", "w");
 kmax = fopen("kmax.out", "w");
 cmax = fopen("cmax.out", "w");
 cmax_t = fopen("cmax_t.out", "w");
 cmax_first = fopen("cmax_first.out", "w");
 cmax_first_t = fopen("cmax_first_t.out", "w");
 totprob = fopen("totprob.out", "w");
 energy = fopen("energy.out", "w");
 times = fopen("times.out", "w");

 // total population k states
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ik_vib; j < Ic_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  fprintf(tkprob, "%-.9g %-.9g\n", t[i], summ);
 }

 // population k states
 for (i = 0; i <= timesteps; i++) {
  for (j = Ik_vib; j < Ic_vib; j++) {
   if (j == Ik_vib) {
    fprintf(kprobs, "%-.9g", t[i]);
   }
   fprintf(kprobs, " %-.9g", pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2));
  }
  fprintf(kprobs, "\n");
 }

 // total population c states
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ic_vib; j < Ib_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  fprintf(tcprob, "%-.9g %-.9g\n", t[i], summ);
 }

 // population c states
 for (i = 0; i <= timesteps; i++) {
  for (j = Ic_vib; j < Ib_vib; j++) {
   if (j == Ic_vib) {
    fprintf(cprobs, "%-.9g", t[i]);
   }
   fprintf(cprobs, " %-.9g", pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2));
  }
  fprintf(cprobs, "\n");
 }

 // total population b states
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ib_vib; j < Il_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  fprintf(tbprob, "%-.9g %-.9g\n", t[i], summ);
 }

 // population b states
 for (i = 0; i <= timesteps; i++) {
  for (j = Ib_vib; j < Il_vib; j++) {
   if (j == Ib_vib) {
    fprintf(bprobs, "%-.9g", t[i]);
   }
   fprintf(bprobs, " %-.9g", pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2));
  }
  fprintf(bprobs, "\n");
 }

 // total population l states
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Il_vib; j < NEQ_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  fprintf(tlprob, "%-.9g %-.9g\n", t[i], summ);
 }

 // population l states
 for (i = 0; i <= timesteps; i++) {
  for (j = Il_vib; j < NEQ_vib; j++) {
   if (j == Il_vib) {
    fprintf(lprobs, "%-.9g", t[i]);
   }
   fprintf(lprobs, " %-.9g", pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2));
  }
  fprintf(lprobs, "\n");
 }

 // max population on k states
 maxx = 0;
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ik_vib; j < Ic_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  if (summ > maxx) {
   maxx = summ;
  }
 }
 fprintf(kmax, "%-.9g\n", maxx);

 // max population on c states and time index
 maxx = 0;
 maxx_t = 0;
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ic_vib; j < Ib_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  if (summ > maxx) {
   maxx = summ;
   maxx_t = i;
  }
 }
 fprintf(cmax, "%-.9g\n", maxx);
 fprintf(cmax_t, "%-.9g\n", t[maxx_t]);

 // first max population on c states and time index
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = Ic_vib; j < Ib_vib; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  if (j == Ic_vib) {
   maxx = summ;
   maxx_t = i;
  }
  else {
   if (summ < maxx) {
    fprintf(cmax_first, "%-.9g\n", maxx);
    fprintf(cmax_first_t, "%-.9g\n", t[maxx_t]);
    break;
   }
   if (summ > maxx) {
    maxx = summ;
    maxx_t = i;
   }
   // if the last data point equals the first maximum
   // (if the population on c is static, basically)
   if ((summ == maxx) && (i == (timesteps - 1))) {
    fprintf(cmax_first, "%-.9g\n", maxx);
    fprintf(cmax_first_t, "%-.9g\n", t[maxx_t]);
   }
  }
  if (summ < maxx) {
   break;
  }
 }

 // total population on all states
 for (i = 0; i <= timesteps; i++) {
  summ = 0;
  for (j = 0; j < dim; j++) {
   summ += pow(psi_t[i*dim + j].re,2) + pow(psi_t[i*dim + j].im,2);
  }
  fprintf(totprob, "%-.9g %-.9g\n", t[i], summ);
 }

 fclose(tkprob);
 fclose(tcprob);
 fclose(tbprob);
 fclose(tlprob);
 fclose(kprobs);
 fclose(cprobs);
 fclose(bprobs);
 fclose(lprobs);
 fclose(kmax);
 fclose(cmax);
 fclose(cmax_t);
 fclose(cmax_first);
 fclose(cmax_first_t);
 fclose(totprob);
 fclose(energy);
 fclose(times);
}

int main (int argc, char * argv[]) {

 // VARIABLES GO HERE//
 int i, j;					// counter!
 int flag;
 realtype * k_pops;				// pointers to arrays of populations
 realtype * l_pops;
 realtype * c_pops;
 realtype * b_pops;
 realtype * ydata;				// pointer to ydata (contains all populations)
 realtype * k_energies;				// pointers to arrays of energies
 realtype * c_energies;
 realtype * b_energies;
 realtype * l_energies;
 int Nk_first;					// first k state initially populated
 int Nk_final;					// final k state initially populated
 realtype bulk_gap;				// bulk band gap
 double valenceBand;				// valence band width
 double temperature;				// system temperature
 double bulkGaussSigma;				// width of initial Gaussian in bulk
 double bulkGaussMu;				// position of initial Gaussian above band edge
 realtype t0 = 0.0;				// initial time
 realtype t = 0;
 realtype tret;					// time returned by the solver
 time_t startRun;				// time at start of log
 time_t endRun;					// time at end of log
 struct tm * currentTime;			// time structure for localtime
#ifdef DEBUG
 FILE * realImaginary;				// file containing real and imaginary parts of the wavefunction
#endif
 FILE * log;					// log file with run times
 realtype * tkprob; 				// total probability in k, l, c, b states at each timestep
 realtype * tlprob;
 realtype * tcprob;
 realtype * tbprob;
 realtype ** vibprob;
 double ** allprob;				// populations in all states at all times
 realtype * times;
 realtype * qd_est;
 realtype * qd_est_diag;
 realtype * energy_expectation;			// expectation value of energy at each timestep
 // END VARIABLES //

 // OPEN LOG FILE; PUT IN START TIME //
 log = fopen("log.out", "w");			// note that this file is closed at the end of the program
 time(&startRun);
 currentTime = localtime(&startRun);
 fprintf(log, "Run started at %s\n", asctime(currentTime));
 
 // read in parameters from parameter bash script

 // ASSIGN VARIABLE DEFAULTS //
 i = 0;
 double summ = 0;			// sum variable
 bool timedepH = 1;			// if H is TD, use CVODE, else diag H and propogate
 bool analytical = 0;			// turn on analytical propagation
 realtype abstol = 1e-10;		// absolute tolerance (for SUNDIALS)
 realtype reltol = 1e-10;		// relative tolerance (for SUNDIALS)
 realtype tout = 10000;			// final time reached by solver in atomic units
 int numsteps = 10000;			// number of time steps
 numOutputSteps = 1000;
 // bulk parameters //
 k_bandedge = 0.0;			// lower band edge of conduction band
 k_bandtop = 0.01;			// upper band edge of bulk conduction band
 bulk_gap = 0.001;			// bulk band gap
 valenceBand = 0.01;			// valence band width
 Nk = 100;				// number of k states
 Nk_first = 1;				// first k state initially populated
 Nk_final = 1;				// final k state initially populated
 bulkGaussSigma = 0.001;		// width of initial Gaussian in bulk
 bulkGaussMu = 0.01;			// position of initial Gaussian above band edge
 // physical parameters //
 temperature = 3e2;			// temperature of the system
 // vibronic parameters //
 N_vib = 1;				// number of vibronic states
 E_vib = 0.001;				// vibrational energy
 gkc = 0.0;				// g factor between k and c states
 gkb = 0.0;				// g factor between k and b states
 gbc = 0.0;				// g factor between b and c states
 gbb = 0.0;				// g factor between b states
 muLK = 1.0;				// transition dipole moment from l to k (energy a.u.)
 pumpFWHM = 1000;
 pumpPeak = 2000;
 pumpFreq = 0.01;
 pumpAmpl = 1.0;
 pumpPhase = 0.0;
 // DONE ASSIGNING VARIABLE DEFAULTS //

 string line;
 string input_param;
 string param_val;
 size_t equals_pos;
 size_t space_pos;

 ifstream bash_in;	// declare input file stream

 bash_in.open("ins/parameters.in", ios::in);	// open file as input stream
 if (bash_in.good() == false) {
  fprintf(stderr, "ERROR [Inputs]: file 'ins/parameters.in' not available for reading\n");
  return -1;
 }

 cout << endl;

 // read first line of input file
 getline (bash_in,line);

 // skip non-parameter lines
 while ( line != "## START INPUT PARAMETERS ##") {
  getline (bash_in,line);
 }

 while ( line != "## END INPUT PARAMETERS ##") {
  // skip comment lines
  if ( line.substr(0,1) == "#" ) {
   getline (bash_in,line);
   continue;
  }
  // find first equals sign
  equals_pos=line.find("=");
  // find first whitespace
  space_pos=(line.find(" ") > line.find("\t") ? line.find("\t") : line.find(" "));
  // parameter name is before equals sign
  input_param = line.substr(0,int(equals_pos));
  // parameter is after equals sign, before space
  param_val = line.substr(int(equals_pos)+1,int(space_pos)-int(equals_pos));
  // extract parameters
#ifdef DEBUG
  cout << "Parameter: " << input_param << endl << "New value: " << atof(param_val.c_str()) << endl;
#endif
  if (input_param == "timedepH") { timedepH = atoi(param_val.c_str()); }
  else if (input_param == "analytical") { analytical = atof(param_val.c_str()); }
  else if (input_param == "abstol") { abstol = atof(param_val.c_str()); }
  else if (input_param == "reltol" ) { reltol = atof(param_val.c_str()); }
  else if (input_param == "tout" ) { tout = atof(param_val.c_str()); }
  else if (input_param == "numsteps" ) { numsteps = atoi(param_val.c_str()); }
  else if (input_param == "numOutputSteps" ) { numOutputSteps = atoi(param_val.c_str()); }
  else if (input_param == "k_bandedge" ) { k_bandedge = atof(param_val.c_str()); }
  else if (input_param == "k_bandtop" ) { k_bandtop = atof(param_val.c_str()); }
  else if (input_param == "bulk_gap" ) { bulk_gap = atof(param_val.c_str()); }
  else if (input_param == "Nk" ) { Nk = atoi(param_val.c_str()); }
  else if (input_param == "Nk_first" ) { Nk_first = atoi(param_val.c_str()); }
  else if (input_param == "Nk_final" ) { Nk_final = atoi(param_val.c_str()); }
  else if (input_param == "valenceBand" ) { valenceBand = atof(param_val.c_str()); }
  else if (input_param == "Nl" ) { Nl = atoi(param_val.c_str()); }
  else if (input_param == "bulkGaussSigma" ) { bulkGaussSigma = atof(param_val.c_str()); }
  else if (input_param == "bulkGaussMu" ) { bulkGaussMu = atof(param_val.c_str()); }
  else if (input_param == "temperature" ) { temperature = atof(param_val.c_str()); }
  else if (input_param == "N_vib" ) { N_vib = atoi(param_val.c_str()); }
  else if (input_param == "E_vib" ) { E_vib = atof(param_val.c_str()); }
  else if (input_param == "gkc" ) { gkc = atof(param_val.c_str()); }
  else if (input_param == "gkb" ) { gkb = atof(param_val.c_str()); }
  else if (input_param == "gbc" ) { gbc = atof(param_val.c_str()); }
  else if (input_param == "gbb" ) { gbb = atof(param_val.c_str()); }
  else if (input_param == "muLK" ) { muLK = atof(param_val.c_str()); }
  else if (input_param == "pumpFWHM" ) { pumpFWHM = atof(param_val.c_str()); }
  else if (input_param == "pumpPeak" ) { pumpPeak = atof(param_val.c_str()); }
  else if (input_param == "pumpFreq" ) { pumpFreq = atof(param_val.c_str()); }
  else if (input_param == "pumpAmpl" ) { pumpAmpl = atof(param_val.c_str()); }
  else if (input_param == "pumpPhase" ) { pumpPhase = atof(param_val.c_str()); }
  else if (input_param == "bulk_FDD" ) { bulk_FDD = atoi(param_val.c_str()); }
  else if (input_param == "bulk_Gauss" ) { bulk_Gauss = atoi(param_val.c_str()); }
  else if (input_param == "bulk_constant" ) { bulk_constant = atoi(param_val.c_str()); }
  else if (input_param == "qd_pops" ) { qd_pops = atoi(param_val.c_str()); }
  else if (input_param == "laser_on" ) { laser_on = atoi(param_val.c_str()); }
  else if (input_param == "parabolicCoupling" ) { parabolicCoupling = atoi(param_val.c_str()); }
  else if (input_param == "scale_bubr" ) { scale_bubr = atoi(param_val.c_str()); }
  else if (input_param == "scale_brqd" ) { scale_brqd = atoi(param_val.c_str()); }
  else if (input_param == "scale_buqd" ) { scale_buqd = atoi(param_val.c_str()); }
  else if (input_param == "scale_laser" ) { scale_laser = atoi(param_val.c_str()); }
  else if (input_param == "bridge_on" ) { bridge_on = atoi(param_val.c_str()); }
  else if (input_param == "random_phase" ) { random_phase = atoi(param_val.c_str()); }
  else if (input_param == "random_seed" ) { random_seed = atoi(param_val.c_str()); }
  else {  }
  getline (bash_in,line);
 }
#ifdef DEBUG
 cout << endl;
 cout << "timedepH is " << timedepH << endl;
 cout << "analytical is " << analytical << endl;
 cout << "abstol is " << abstol << endl;
 cout << "reltol is " << reltol << endl;
 cout << "tout is " << tout << endl;
 cout << "numsteps is " << numsteps << endl;
 cout << "numOutputSteps is " << numOutputSteps << endl;
 cout << "k_bandedge is " << k_bandedge << endl;
 cout << "k_bandtop is " << k_bandtop << endl;
 cout << "bulk_gap is " << bulk_gap << endl;
 cout << "Nk is " << Nk << endl;
 cout << "Nk_first is " << Nk_first << endl;
 cout << "Nk_final is " << Nk_final << endl;
 cout << "valenceBand is " << valenceBand << endl;
 cout << "Nl is " << Nl << endl;
 cout << "bulkGaussSigma is " << bulkGaussSigma << endl;
 cout << "bulkGaussMu is " << bulkGaussMu << endl;
 cout << "temperature is " << temperature << endl;
 cout << "N_vib is " << N_vib << endl;
 cout << "E_vib is " << E_vib << endl;
 cout << "gkc is " << gkc << endl;
 cout << "gkb is " << gkb << endl;
 cout << "gbc is " << gbc << endl;
 cout << "gbb is " << gbb << endl;
 cout << "muLK is " << muLK << endl;
 cout << "pumpFWHM is " << pumpFWHM << endl;
 cout << "pumpPeak is " << pumpPeak << endl;
 cout << "pumpFreq is " << pumpFreq << endl;
 cout << "pumpAmpl is " << pumpAmpl << endl;
 cout << "pumpPhase is " << pumpPhase << endl;
 cout << "bulk_FDD is " << bulk_FDD << endl;
 cout << "bulk_Gauss is " << bulk_Gauss << endl;
 cout << "bulk_constant is " << bulk_constant << endl;
 cout << "qd_pops is " << qd_pops << endl;
 cout << "laser_on is " << laser_on << endl;
 cout << "parabolicCoupling is " << parabolicCoupling << endl;
 cout << "scale_bubr is " << scale_bubr << endl;
 cout << "scale_brqd is " << scale_brqd << endl;
 cout << "scale_buqd is " << scale_buqd << endl;
 cout << "scale_laser is " << scale_laser << endl;
 cout << "bridge_on is " << bridge_on << endl;
 cout << "random_phase is " << random_phase << endl;
 cout << "random_seed is " << random_seed << endl;
#endif

 // make a note about the laser intensity.
 fprintf(log,"The laser intensity is %.5e W/cm^2.\n\n",pow(pumpAmpl,2)*3.5094452e16);

 // Error checking
 if ((bulk_FDD && qd_pops) || (bulk_constant && qd_pops) || (bulk_Gauss && qd_pops)) {
  cerr << "\nWARNING: population starting both in bulk and QD.\n";
 }
 if (Nk_first > Nk || Nk_first < 1) {
  fprintf(stderr, "ERROR [Inputs]: Nk_first greater than Nk or less than 1.\n");
  return -1;
 }
 if (Nk_final > Nk || Nk_final < 1) {
  fprintf(stderr, "ERROR [Inputs]: Nk_final greater than Nk or less than 1.\n");
  return -1;
 }
 if (Nk_final < Nk_first) {
  fprintf(stderr, "ERROR [Inputs]: Nk_final is less than Nk_first.\n");
  return -1;
 }
 if (Nl < 1) {
  fprintf(stderr, "ERROR [Inputs]: Nl less than 1.\n");
  return -1;
 }
 if ((bulk_FDD && bulk_constant) || (bulk_FDD && bulk_Gauss) || (bulk_constant && bulk_Gauss)) {
  cerr << "\nERROR: two different switches are on for bulk starting conditions.\n";
  return -1;
 }
 if (random_seed < -1) {
  cerr << "\nERROR: random_phase must be -1 or greater.\n";
  return -1;
 }

 cout << endl;

 bash_in.close();

 // DONE ASSIGNING VARIABLES FROM RUN SCRIPT //

 // READ DATA FROM INPUTS //
 Nc = Number_of_values("ins/c_energies.in");
 Nb = Number_of_values("ins/b_energies.in");
 k_pops = new realtype [Nk];
 c_pops = new realtype [Nc];
 b_pops = new realtype [Nb];
 l_pops = new realtype [Nl];
 k_energies = new realtype [Nk];
 c_energies = new realtype [Nc];
 b_energies = new realtype [Nb];
 l_energies = new realtype [Nl];
 if (Number_of_values("ins/c_pops.in") != Nc)
  fprintf(stderr, "ERROR [Inputs]: c_pops and c_energies not the same length.\n");
 Read_array_from_file(c_energies, "ins/c_energies.in", Nc);
 if (bridge_on) {
  if (bridge_on && (Nb < 1)) {
   cerr << "\nERROR: bridge_on but no bridge states.  The file b_energies.in is probably empty.\n";
   return -1;
  }
  Vbridge = new realtype [Nb+1];
  Read_array_from_file(b_energies, "ins/b_energies.in", Nb);
  Read_array_from_file(Vbridge, "ins/Vbridge.in", Nb + 1);
 }
 else {
  Nb = 0;
  Vnobridge = new realtype [1];
  Read_array_from_file(Vnobridge, "ins/Vnobridge.in", 1);
 }
 // DONE READING //
#ifdef DEBUG
 cout << "\nDone reading things from inputs.\n";
#endif

 // PREPROCESS DATA FROM INPUTS //
 NEQ = Nk+Nc+Nb+Nl;				// total number of equations set
 NEQ_vib = NEQ*N_vib;
#ifdef DEBUG
 cout << "\nTotal number of states: " << NEQ << endl;
 cout << Nk << " bulk, " << Nc << " QD, " << Nb << " bridge, " << Nl << " bulk VB.\n";
#endif
 tkprob = new realtype [numOutputSteps+1];	// total population on k, b, c at each timestep
 tcprob = new realtype [numOutputSteps+1];
 tbprob = new realtype [numOutputSteps+1];
 tlprob = new realtype [numOutputSteps+1];
 vibprob = new realtype * [numOutputSteps+1];
 for (i = 0; i < numOutputSteps+1; i++)
  vibprob[i] = new realtype [N_vib];
 allprob = new double * [numOutputSteps+1];
 for (i = 0; i < numOutputSteps+1; i++)
  allprob[i] = new double [NEQ];
 times = new realtype [numOutputSteps+1];
 qd_est = new realtype [numOutputSteps+1];
 qd_est_diag = new realtype [numOutputSteps+1];
 energy_expectation = new realtype [numOutputSteps+1];	// expectation value of energy; for sanity checking
 Ik = 0;					// set index start positions for each type of state
 Ic = Nk;
 Ib = Ic+Nc;
 Il = Ib+Nb;
 Ik_vib = 0;
 Ic_vib = Nk*N_vib;
 Ib_vib = Ic_vib + Nc*N_vib;
 Il_vib = Ib_vib + Nb*N_vib;
 // assign bulk conduction band energies
 Build_continuum(k_energies, Nk, k_bandedge, k_bandtop);
 // assign bulk valence band energies
 Build_continuum(l_energies, Nl, k_bandedge - valenceBand - bulk_gap, k_bandedge - bulk_gap);
 // assign populations
 Initialize_array(b_pops, Nb, 0.0);		// populate b states
 if (bulk_FDD) {
  Build_k_pops(k_pops, k_energies, k_bandedge, temperature);   // populate k states with FDD
  Initialize_array(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  Initialize_array(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (bulk_constant) {
#ifdef DEBUG
  cout << "\ninitializing k_pops\n";
#endif
  Initialize_array(k_pops, Nk, 0.0);
#ifdef DEBUG
  cout << "\ninitializing k_pops with constant population in states\n";
#endif
  Initialize_array(k_pops+Nk_first-1, Nk_final-Nk_first+1, 1.0);
#ifdef DEBUG
  cout << "\nThis is k_pops:\n";
  for (i = 0; i < Nk; i++) {
   cout << k_pops[i] << endl;
  }
  cout << "\n";
#endif
  Initialize_array(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  Initialize_array(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (bulk_Gauss) {
  Build_k_pops_Gaussian(k_pops, k_energies, k_bandedge,
                        bulkGaussSigma, bulkGaussMu);   // populate k states with FDD
  Initialize_array(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  Initialize_array(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (qd_pops) {
  Read_array_from_file(c_pops, "ins/c_pops.in", Nc);	// QD populations from file
  Initialize_array(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  Initialize_array(k_pops, Nk, 0.0);             // populate k states (all zero to start off)
 }
 else {
  Initialize_array(k_pops, Nk, 0.0);             // populate k states (all zero to start off)
  Initialize_array(l_pops, Nl, 1.0);		// populate l states (all populated to start off)
  Initialize_array(c_pops, Nc, 0.0);		// QD states empty to start
 }
 if (Nb == 0) {					// assign Franck-Condon factors
  FCkc = new realtype * [N_vib];
  for (i = 0; i < N_vib; i++)
   FCkc[i] = new realtype [N_vib];
  Build_Franck_Condon_factors(FCkc, gkc, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCkc:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7e ", FCkc[i][j]);
    cout << endl;
   }
#endif
 }
 else if (bridge_on) {
  if (Nb > 1) {
   FCbb = new realtype * [N_vib];
   for (i = 0; i < N_vib; i++)
    FCbb[i] = new realtype [N_vib];
   Build_Franck_Condon_factors(FCbb, gbb, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCbb:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7g ", FCbb[i][j]);
    cout << endl;
   }
#endif
  }
  FCkb = new realtype * [N_vib];
  for (i = 0; i < N_vib; i++)
   FCkb[i] = new realtype [N_vib];
  Build_Franck_Condon_factors(FCkb, gkb, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCkb:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7g ", FCkb[i][j]);
    cout << endl;
   }
#endif
  FCbc = new realtype * [N_vib];
  for (i = 0; i < N_vib; i++)
   FCbc[i] = new realtype [N_vib];
  Build_Franck_Condon_factors(FCbc, gbc, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCbc:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7g ", FCbc[i][j]);
    cout << endl;
   }
#endif
 }

 ydata = new realtype [2*NEQ_vib];			// assemble ydata
 Initialize_array(ydata, 2*NEQ_vib, 0.0);
 for (i = 0; i < Nk; i++)
  ydata[Ik_vib + i*N_vib] = k_pops[i];
 for (i = 0; i < Nc; i++)
  ydata[Ic_vib + i*N_vib] = c_pops[i];
 for (i = 0; i < Nb; i++)
  ydata[Ib_vib + i*N_vib] = b_pops[i];
 for (i = 0; i < Nl; i++)
  ydata[Il_vib + i*N_vib] = l_pops[i];

#ifdef DEBUG
 cout << "\n to start, y data is:\n";
 for (i = 0; i < NEQ_vib; i++) {
  cout << ydata[i] << "\n";
 }
#endif

 outputYData(ydata, NEQ_vib);

 // If random_phase is on, give all coefficients a random phase
 if (random_phase) {
  float phi;
  // set the seed
  if (random_seed == -1) { srand(time(NULL)); }
  else { srand(random_seed); }
  for (i = 0; i < NEQ_vib; i++) {
   phi = 2*3.1415926535*(float)rand()/(float)RAND_MAX;
   ydata[i] = ydata[i]*cos(phi);
   ydata[i + NEQ_vib] = ydata[i + NEQ_vib]*sin(phi);
  }
 }
//these lines are a test
 //Initialize_array(ydata, 2*NEQ_vib, 0.0000001);
 // for (i = 0; i < NEQ_vib; i += 2) {
  // ydata[i] = 0.0000001;
 // }
 // Activate the following line to have the electron start on the first bridge
 // ydata[Ib_vib] = 1.0;
#ifdef DEBUG
 cout << endl;
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Re[k(" << i << "," << j << ")] = " << ydata[Ik_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Re[c(" << i << "," << j << ")] = " << ydata[Ic_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Re[b(" << i << "," << j << ")] = " << ydata[Ib_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nl; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Re[l(" << i << "," << j << ")] = " << ydata[Il_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[k(" << i << "," << j << ")] = " << ydata[Ik_vib + i*N_vib + j + NEQ_vib] << endl;
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[c(" << i << "," << j << ")] = " << ydata[Ic_vib + i*N_vib + j + NEQ_vib] << endl;
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[b(" << i << "," << j << ")] = " << ydata[Ib_vib + i*N_vib + j + NEQ_vib] << endl;
 for (i = 0; i < Nl; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[l(" << i << "," << j << ")] = " << ydata[Il_vib + i*N_vib + j + NEQ_vib] << endl;
 cout << endl;
 summ = 0;
 for (i = 0; i < 2*NEQ_vib; i++) {
  summ += pow(ydata[i],2);
 }
 cout << "\nTotal population is " << summ << "\n\n";
#endif
 energy = new realtype [NEQ_vib];			// assemble energy array
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++)
   energy[Ik_vib + i*N_vib + j] = k_energies[i] + E_vib*j;
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++)
   energy[Ic_vib + i*N_vib + j] = c_energies[i] + E_vib*j;
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
   energy[Ib_vib + i*N_vib + j] = b_energies[i] + E_vib*j;
 for (i = 0; i < Nl; i++)
  for (j = 0; j < N_vib; j++)
   energy[Il_vib + i*N_vib + j] = l_energies[i] + E_vib*j;
 user_data = energy;
#ifdef DEBUG
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++)
  cout << "energy[k(" << i << "," << j << ")] = " << energy[Ik_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++)
  cout << "energy[c(" << i << "," << j << ")] = " << energy[Ic_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
  cout << "energy[b(" << i << "," << j << ")] = " << energy[Ib_vib + i*N_vib + j] << endl;
 for (i = 0; i < Nl; i++)
  for (j = 0; j < N_vib; j++)
  cout << "energy[l(" << i << "," << j << ")] = " << energy[Il_vib + i*N_vib + j] << endl;
 cout << endl;
#endif

 // assign coupling constants
 V = new realtype * [NEQ];
 for (i = 0; i < NEQ; i++)
  V[i] = new realtype [NEQ];
 buildCoupling(V, NEQ, k_bandedge, k_bandtop, energy);

 // make a note in the log about system timescales
 double tau = 0;		// fundamental system timescale
 if (Nk == 1) {
  fprintf(log, "\nThe timescale (tau) is undefined (Nk == 1).\n");
 }
 else {
  if (bridge_on) {
   if (scale_bubr) {
    tau = 1.0/(2*Vbridge[0]*M_PI);
   }
   else {
    tau = ((k_bandtop - k_bandedge)/(Nk - 1))/(2*pow(Vbridge[0],2)*M_PI);
   }
  }
  else {
   if (scale_buqd) {
    tau = 1.0/(2*Vnobridge[0]*M_PI);
   }
   else {
    tau = ((k_bandtop - k_bandedge)/(Nk - 1))/(2*pow(Vnobridge[0],2)*M_PI);
   }
  }
  fprintf(log, "\nThe timescale (tau) is %.9e a.u.\n", tau);
 }

 // DONE PREPROCESSING //

 // Creates N_Vector y with initial populations which will be used by CVode//
 y = N_VMake_Serial(2*NEQ_vib, ydata);
 yout = N_VClone(y);

 // print t = 0 information //
 Normalize_NV(y, 1.00);			// normalizes all populations to 1; this is for one electron
 summ = 0;
 for (i = 0; i < 2*NEQ_vib; i++) {
  summ += pow(NV_Ith_S(y, i),2);
 }
#ifdef DEBUG
  cout << "\nAfter normalization, total population is " << summ << "\n\n";
#endif
 if ( summ == 0.0 ) {
  cerr << "\nFATAL ERROR [populations]: total population is 0!\n";
  return -1;
 }
 if ( fabs(summ-1.0) > 1e-12 ) {
  cerr << "\nWARNING [populations]: total population is not 1, it is " << summ << "!\n";
 }
#ifdef DEBUG
 realImaginary = fopen("real_imaginary.out", "w");
#endif
 Output_checkpoint(
#ifdef DEBUG
   realImaginary, 
#endif
   allprob, y, t0, tkprob, tlprob, tcprob, tbprob, vibprob, times, qd_est,
   qd_est_diag, energy_expectation, 0, energy, k_bandedge, k_bandtop, k_pops);

 if (analytical) {
  // Compute the analytical population on the c states
  Analytical_c(tout, numOutputSteps, energy, k_bandedge, k_bandtop, y);
 }

 // create CVode object //
 cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);	// this is a stiff problem, I guess?
 flag = CVodeSetUserData(cvode_mem, (void *) user_data);	// now stuff in energy is available to CVode

 // initialize CVode solver //
 flag = CVodeInit(cvode_mem, &f, t0, y);

 // specify integration tolerances //
 flag = CVodeSStolerances(cvode_mem, reltol, abstol);

 // attach linear solver module //
 flag = CVDense(cvode_mem, 2*NEQ_vib);
#ifdef DEBUG_SAI
 // match up the timesteps with Sai's
 /*flag = CVodeSetInitStep(cvode_mem, (tout - t0)/((double) numsteps));
 if (flag == CV_MEM_NULL)
  fprintf(stderr, "ERROR [CVodeSetInitStep]: cvode_mem pointer is null");*/
 cout << (tout - t0)/((double) numsteps) << endl;
 flag = CVodeSetMinStep(cvode_mem, (tout - t0)/((double) numsteps));
 if (flag == CV_MEM_NULL)
  fprintf(stderr, "ERROR [CVodeSetMinStep]: cvode_mem pointer is null");
 if (flag == CV_ILL_INPUT)
  fprintf(stderr, "ERROR [CVodeSetMinStep]: hmin is nonpositive or it exceeds the maximum allowable step size.");
 flag = CVodeSetMaxStep(cvode_mem, (tout - t0)/((double) numsteps));
 if (flag == CV_MEM_NULL)
  fprintf(stderr, "ERROR [CVodeSetMaxStep]: cvode_mem pointer is null");
 if (flag == CV_ILL_INPUT)
  fprintf(stderr, "ERROR [CVodeSetMaxStep]: hmin is nonpositive or it exceeds the maximum allowable step size.");
 // specify integration tolerances; don't much care about errors here
 flag = CVodeSStolerances(cvode_mem, 1.0e-1, 1.0e-1);
#endif

 // advance the solution in time! //
 // use CVODE for time-dependent H
 if (timedepH) {
  for (i = 1; i <= numsteps; ++i) {
   t = (tout*((double) i)/((double) numsteps));
   flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
   cout << endl << "CVode flag at step " << i << ": " << flag << endl;
#endif
   if (i % (numsteps/numOutputSteps) == 0) {
    fprintf(stderr, "\r%-.2lf percent done", ((double)i/((double)numsteps))*100);
    Output_checkpoint(
#ifdef DEBUG
      realImaginary, 
#endif
      allprob, yout, t, tkprob, tlprob, tcprob, tbprob, vibprob, times, qd_est,
      qd_est_diag, energy_expectation, (i*numOutputSteps/numsteps), energy,
      k_bandedge, k_bandtop, k_pops);
   }
  }

 // compute final outputs //
 Compute_final_outputs(allprob, times, tkprob,
   tlprob, tcprob, tbprob, vibprob, energy,
   energy_expectation, numOutputSteps, qd_est, qd_est_diag);
 }
 // use LAPACK for time-independent H
 else {
#ifdef DEBUG
  fprintf(stderr, "calculating output times.\n");
#endif
  // calculate output times
  for (i = 0; i <= numOutputSteps; i++) {
   times[i] = tout*((double)i/((double)numOutputSteps));
  }
  // build Hamiltonian
#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif
  realtype * H = new realtype [NEQ_vib*NEQ_vib];
  buildHamiltonian(H, energy, V, NEQ_vib, N_vib, FCkb, FCbb, FCbc);
  printSquareMatrix(H, NEQ_vib, "ham.out");
  // declare LAPACK variables
  char JOBZ;            // 'N' to just compute evals; 'V' to compute evals and evecs
  char UPLO;            // 'U' to store upper diagonal of matrix
  int N;                // order of matrix
  int LDA;              // leading dimension of matrix
  double * W;           // output array of evals
  double * WORK;	// working memory to use during diagonalization
  int LWORK;            // length of WORK; LWORK >= (MAX(1,LWORK))
  int INFO;             // state flag: (0 success, -i ith argument bad, >0 fail)
  // assign LAPACK variables
  JOBZ = 'V';
  UPLO = 'U';
  N = NEQ_vib;
  LDA = N;
  W = new double [N];
  LWORK = -1;
  WORK = new double [1];
  INFO = 0;
  // dry run to compute LWORK
  // using same arguments as actual calculation
#ifdef DEBUG
  fprintf(stderr, "Diagonalizing Hamiltonian.\n");
#endif
  dsyev_(&JOBZ, &UPLO, &N, H, &LDA, W, WORK, &LWORK, &INFO);
  // assign LWORK from first element of WORK
  LWORK = WORK[0];
  // reallocate WORK
  delete [] WORK;
  WORK = new double [LWORK];
  // diagonalize the guy!
  dsyev_(&JOBZ, &UPLO, &N, H, &LDA, W, WORK, &LWORK, &INFO);
  // print eigenvalues and eigenvectors
  printVector(W, N, "evals.out");
  printSquareMatrix(H, N, "evecs.out");
  // make a complex array to represent the starting psi (site basis)
  complex16 * psi_S = new complex16 [NEQ_vib];
  for (i = 0; i < NEQ_vib; i++) {
   psi_S[i].re = NV_Ith_S(y, i);
   psi_S[i].im = NV_Ith_S(y, i+NEQ_vib);
  }
  // make a complex array to represent the starting psi (eigenstate basis)
  complex16 * psi_E = new complex16 [NEQ_vib];
  // project the starting wavefunction onto the eigenstate basis
  projectSiteToState(psi_S, NEQ_vib, H, psi_E);
  // print the starting wavefunction in the two bases
  printCVector(psi_S, NEQ_vib, "psi_start_s.out");
  printCVector(psi_E, NEQ_vib, "psi_start_e.out");
  printPsiSquare(psi_E, W, NEQ_vib, "psi2_start_e.out");
  // make arrays to represent the wavefunction in time
  complex16 * psi_S_t = new complex16 [NEQ_vib*(numOutputSteps+1)];
  complex16 * psi_E_t = new complex16 [NEQ_vib*(numOutputSteps+1)];
  // propagate the wavefunction in time
  propagatePsi(psi_E, psi_E_t, NEQ_vib, W, numOutputSteps, tout);
  // print out the propagated wavefunction (eigenstate basis)
  printCVectorTime(psi_E_t, NEQ_vib, (numOutputSteps+1), "psi_e_t.out");
  // project back onto the site basis
  projectStateToSite(psi_E_t, NEQ_vib, H, psi_S_t, numOutputSteps);
  // print out the propagated wavefunction (site basis)
  printCVectorTime(psi_S_t, NEQ_vib, (numOutputSteps+1), "psi_s_t.out");
  makeOutputsTI(psi_S_t, NEQ_vib, times, numOutputSteps);
  // write out projections of subsystems
  projectSubsystems(H, W, NEQ_vib);
  delete [] H;
  delete [] W;
  delete [] WORK;
  delete [] psi_S;
  delete [] psi_E;
  delete [] psi_S_t;
  delete [] psi_E_t;
 }

 // compute time-independent outputs
 FILE * energyFile = fopen("energy.out", "w");
 for (i = 0; i < NEQ_vib; i++) {
  fprintf(energyFile, "%-.9e\n", energy[i]);
 }
 fclose(energyFile);
#ifdef DEBUG
 fclose(realImaginary);
#endif
 
 // finalize log file //
 time(&endRun);
 currentTime = localtime(&endRun);
 fprintf(log, "Final status of 'flag' variable: %d\n\n", flag);
 fprintf(log, "Run ended at %s\n", asctime(currentTime));
 fprintf(log, "Run took %.3g seconds.\n", difftime(endRun, startRun));
 fclose(log);					// note that the log file is opened after variable declaration
 printf("\nRun took %.3g seconds.\n", difftime(endRun, startRun));

 // deallocate memory for N_Vectors //
 N_VDestroy_Serial(y);
 N_VDestroy_Serial(yout);

 // free solver memory //
 CVodeFree(&cvode_mem);

 // delete all these guys
 delete [] tkprob;
 delete [] tlprob;
 delete [] tcprob;
 delete [] tbprob;
 delete [] vibprob;
 delete [] k_pops;
 delete [] c_pops;
 delete [] b_pops;
 delete [] energy;
 delete [] V;
 delete [] Vbridge;
 delete [] Vnobridge;
 delete [] k_energies;
 delete [] c_energies;
 delete [] b_energies;
 delete [] l_energies;
 delete [] ydata;
 delete [] FCkc;
 delete [] FCkb;
 delete [] FCbc;
 delete [] FCbb;
 fprintf(stderr, "\nwhoo\n");
 return 0;
}