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

/* DEBUG compiler flag: turn this on to generate basic debug outputs.         */
//#define DEBUG
/* DANGER! Only turn on DEBUGf for small test runs, otherwise output is       */
/* enormous (many GB).  This flag turns on debug output within the f          */
/* function.                                                                  */
//#define DEBUGf
/* This flag is debuggery related to checking against Sai's code.             */
//#define DEBUG_SAI

using namespace std;

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
 int bulk_FDD = 0;			// switches for starting conditions
 int bulk_Gauss = 0;
 int bulk_constant = 0;
 int qd_pops = 0;
 int laser_on = 0;
 int scale_bubr = 0;
 int scale_brqd = 0;
 int scale_buqd = 0;
 int scale_laser = 0;
 int bridge_on = 0;
 int random_phase = 0;
 int random_seed = 0;
// END GLOBAL VARIABLES
#ifdef DEBUG_SAI
 double last_t = -1.0;				// keeps track of last time for which debuggery was printed
#endif

// returns the number of numbers in a file.  This way, it doesn't matter if
// they are one per line or multiple per line.
int Number_of_values (const char * nameOfFile) {
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

// reads in the values from file; returns an array the length of the number of 
// numbers in the file
void Read_array_from_file (realtype * array, const char * nameOfFile, int numberOfValues) {

 FILE * input;
 int i = 0;

 input = fopen(nameOfFile,"r");

 if (input != NULL) {
  while (fscanf(input, "%lf", &array[i]) != EOF) {
   i++;
  }
 }
 else {
  fprintf(stderr, "WARNING [Read_array_from_file]: file %s does not exist.\n", nameOfFile);
 }

 fclose(input);
}

// Returns an array of length n with all values set to initializeValue. //
void Initialize_array(realtype * array, int n, realtype initializeValue) {

 int i;

 for (i = 0; i < n; i++) {
  array[i] = initializeValue;
 }
}


void Build_continuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
 
 int i;

 Energies[0] = BandEdge;	// the bottom of the conduction band is set
 
 // loop over the remaining states.  This way the top of the band will be at BandTop
 for (i = 1; i < numberOfStates; i++) {
  Energies[i] = Energies[i-1] + (BandTop-BandEdge)/(numberOfStates-1);
 }
}


void Build_k_pops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp) {

 int i;

 for (i = 0; i < Nk; i++) {
  kPops[i] = sqrt(1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.01)*3.185e5/(temp))));	// dunno where the actual Fermi level is
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

// assign coupling constants to global array V //
void Build_v (realtype ** vArray, int dim, realtype kBandEdge, realtype kBandTop) {
 
 int i, j;					// counters

 if ((scale_buqd) && (Nk > 1))
  Vnobridge[0] = sqrt(Vnobridge[0]*(kBandTop-kBandEdge)/(Nk-1));

 for (i = 0; i < dim; i++)			// initialize
  for (j = 0; j < dim; j++)
   vArray[i][j] = 0.0;


 if (bridge_on) {				// bridge
  // coupling between k and b1
  if ((scale_bubr) && (Nk > 1)) {
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = sqrt(Vbridge[0]*(kBandTop-kBandEdge)/(Nk-1));
    vArray[Ib][Ik+i] = sqrt(Vbridge[0]*(kBandTop-kBandEdge)/(Nk-1));
   }
  }
  else {
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = Vbridge[0];
    vArray[Ib][Ik+i] = Vbridge[0];
   }
  }
  // coupling between bN and c
  if ((scale_brqd) && (Nc > 1)) {
   for (i = 0; i < Nc; i++) {
    vArray[Ic+i][Ib+Nb-1] = Vbridge[Nb]/sqrt(Nc-1);
    vArray[Ib+Nb-1][Ic+i] = Vbridge[Nb]/sqrt(Nc-1);
   }
  }
  else {
   for (i = 0; i < Nc; i++) {
    vArray[Ic+i][Ib+Nb-1] = Vbridge[Nb];
    vArray[Ib+Nb-1][Ic+i] = Vbridge[Nb];
   }
  }
  // coupling between bridge states
  for (i = 0; i < Nb - 1; i++) {
   vArray[Ib+i][Ib+i+1] = Vbridge[i+1];
   vArray[Ib+i+1][Ib+i] = Vbridge[i+1];
  }
 }
 else {					// no bridge
  // this is goofy, I should probably not have the sqrt(27.211) scaling factor
  for (i = 0; i < Nk; i++) {
   for (j = 0; j < Nc; j++) {
    vArray[Ik+i][Ic+j] = Vnobridge[0];
    vArray[Ic+j][Ik+i] = Vnobridge[0];
   }
  }
  for (i = 0; i < Nk; i++) {
   for (j = 0; j < Nc; j++) {
    vArray[Ik+i][Ic+j] = Vnobridge[0]/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);
    vArray[Ic+j][Ik+i] = Vnobridge[0]/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);
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
 double sigma = pumpFWHM/2.35482005;
 return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}


// gives f(y,t) //
int f(realtype t, N_Vector y, N_Vector ydot, void * data) {
 
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


// normalizes an N_Vector so that the populations of all states are normalized to value /total/ //
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


int Output_checkpoint(
#ifdef DEBUG
  FILE * realImaginary, 
#endif
  double ** allprobs, N_Vector outputData, realtype time,
  realtype * totK, realtype * totL, realtype * totC, realtype * totB, realtype ** vibProb, realtype * times,
  realtype * qd_est, realtype * qd_est_diag, realtype * energy_expectation, int index, realtype * energies,
  double kBandEdge, double kBandTop, double * k_pops) {

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
 fprintf(realImaginary, "%-.7lf ", time);
 for (i = 0; i < NEQ_vib; i++) {
  fprintf(realImaginary, "%-7lf %-7lf", NV_Ith_S(outputData, i), NV_Ith_S(outputData, i+NEQ_vib));
  if (i == (NEQ_vib - 2)) {
   fprintf(realImaginary, " ");
  }
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

 int i;
 realtype riemann = 0;

 for (i = 0; i < num-1; i++)
  riemann += (values[i+1] + values[i])*(time[i+1]-time[i])/2;

 return riemann;
}


realtype Find_array_maximum (realtype * inputArray, int num) {

 int i;
 realtype currentMax = inputArray[0];

 for (i = 1; i < num; i++)
  if (inputArray[i] > currentMax)
   currentMax = inputArray[i];

 return currentMax;
}


realtype Find_first_array_maximum (realtype * inputArray, int num) {

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
 //
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
 int Nk_init;					// number of k states initially populated
 realtype bulk_gap;				// bulk band gap
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
 realtype abstol = 1e-10;		// absolute tolerance (for SUNDIALS)
 realtype reltol = 1e-10;		// relative tolerance (for SUNDIALS)
 realtype tout = 10000;			// final time reached by solver in atomic units
 int numsteps = 10000;			// number of time steps
 numOutputSteps = 1000;
 // bulk parameters //
 k_bandedge = 0.0;			// lower band edge of conduction band
 k_bandtop = 0.01;			// upper band edge of bulk conduction band
 bulk_gap = 0.001;			// bulk band gap
 Nk = 100;				// number of k states
 Nk_init = 1;				// number of k states initially populated
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

 bash_in.open("ins/parameters.sh", ios::in);	// open file as input stream
 if (bash_in.good() == false) {
  fprintf(stderr, "ERROR [Inputs]: file 'total_dynamix' not available for reading\n");
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
  if (input_param == "abstol") { abstol = atof(param_val.c_str()); }
  else if (input_param == "reltol" ) { reltol = atof(param_val.c_str()); }
  else if (input_param == "tout" ) { tout = atof(param_val.c_str()); }
  else if (input_param == "numsteps" ) { numsteps = atoi(param_val.c_str()); }
  else if (input_param == "numOutputSteps" ) { numOutputSteps = atoi(param_val.c_str()); }
  else if (input_param == "k_bandedge" ) { k_bandedge = atof(param_val.c_str()); }
  else if (input_param == "k_bandtop" ) { k_bandtop = atof(param_val.c_str()); }
  else if (input_param == "bulk_gap" ) { bulk_gap = atof(param_val.c_str()); }
  else if (input_param == "Nk" ) { Nk = atoi(param_val.c_str()); }
  else if (input_param == "Nk_init" ) { Nk_init = atoi(param_val.c_str()); }
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
 cout << "abstol is " << abstol << endl;
 cout << "reltol is " << reltol << endl;
 cout << "tout is " << tout << endl;
 cout << "numsteps is " << numsteps << endl;
 cout << "numOutputSteps is " << numOutputSteps << endl;
 cout << "k_bandedge is " << k_bandedge << endl;
 cout << "k_bandtop is " << k_bandtop << endl;
 cout << "bulk_gap is " << bulk_gap << endl;
 cout << "Nk is " << Nk << endl;
 cout << "Nk_init is " << Nk_init << endl;
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
 if (Nk_init > Nk || Nk_init < 0) {
  fprintf(stderr, "ERROR [Inputs]: Nk_init greater than Nk or less than 0.\n");
  return -1;
 }
 if (bulk_FDD != 0 && bulk_FDD != 1) {
  cerr << "\nERROR: bulk_FDD switch is not 0 or 1.\n";
  return -1;
 }
 if (bulk_Gauss != 0 && bulk_Gauss != 1) {
  cerr << "\nERROR: bulk_Gauss switch is not 0 or 1.\n";
  return -1;
 }
 if (bulk_constant != 0 && bulk_constant != 1) {
  cerr << "\nERROR: bulk_constant switch is not 0 or 1.\n";
  return -1;
 }
 if (qd_pops != 0 && qd_pops != 1) {
  cerr << "\nERROR: qd_pops switch is not 0 or 1.\n";
  return -1;
 }
 if (laser_on != 0 && laser_on != 1) {
  cerr << "\nERROR: laser_on switch is not 0 or 1.\n";
  return -1;
 }
 if (scale_bubr != 0 && scale_bubr != 1) {
  cerr << "\nERROR: scale_bubr switch is not 0 or 1.\n";
  return -1;
 }
 if (scale_brqd != 0 && scale_brqd != 1) {
  cerr << "\nERROR: scale_brqd switch is not 0 or 1.\n";
  return -1;
 }
 if (scale_buqd != 0 && scale_buqd != 1) {
  cerr << "\nERROR: scale_buqd switch is not 0 or 1.\n";
  return -1;
 }
 if (scale_laser != 0 && scale_laser != 1) {
  cerr << "\nERROR: scale_laser switch is not 0 or 1.\n";
  return -1;
 }
 if (bridge_on != 0 && bridge_on != 1) {
  cerr << "\nERROR: bridge_on switch is not 0 or 1.\n";
  return -1;
 }
 if ((bulk_FDD && bulk_constant) || (bulk_FDD && bulk_Gauss) || (bulk_constant && bulk_Gauss)) {
  cerr << "\nERROR: two different switches are on for bulk starting conditions.\n";
  return -1;
 }
 if (random_phase != 0 && random_phase != 1) {
  cerr << "\nERROR: random_phase switch is not 0 or 1.\n";
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
  fprintf(stderr, "ERROR [Inputs]: c_pops and c_energies not the same length.");
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
 Nl = 1;
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
 // assign bulk energies
 Build_continuum(k_energies, Nk, k_bandedge, k_bandtop);	// create bulk conduction quasicontinuum
 l_energies[0]=k_bandedge-bulk_gap;             // assign l energy
 // assign populations
 Initialize_array(b_pops, Nb, 0.0);		// populate b states
 if (bulk_FDD) {
  Build_k_pops(k_pops, k_energies, k_bandedge, temperature);   // populate k states with FDD
  Initialize_array(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  Initialize_array(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (bulk_constant) {
  Initialize_array(k_pops, Nk, 0.0);
  Initialize_array(k_pops, Nk_init, 1.0);
  //k_pops[(Nk/2)+(Nk%2)-1] = 1.0;		// populate just the middle state in the bulk
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
 V = new realtype * [NEQ];
 for (i = 0; i < NEQ; i++)
  V[i] = new realtype [NEQ];
 Build_v(V, NEQ, k_bandedge, k_bandtop);	// assign coupling constants
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

 // If random_phase is on, give all coefficients a random phase
 if (random_phase) {
  float phi;
  // set the seed
  if (random_seed == -1) { srand(time(NULL)); }
  else { srand(random_seed); }
  for (i = 0; i < NEQ_vib; i++) {
   phi = (float)rand()/(float)RAND_MAX;
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
 for (i = 0; i < NEQ_vib; i++) {
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
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
  cout << "energy[l(" << i << "," << j << ")] = " << energy[Il_vib + i*N_vib + j] << endl;
 cout << endl;
#endif
 // DONE PREPROCESSING //

 // Creates N_Vector y with initial populations which will be used by CVode//
 y = N_VMake_Serial(2*NEQ_vib, ydata);
 yout = N_VClone(y);

 // print t = 0 information //
 Normalize_NV(y, 1.00);			// normalizes all populations to 1; this is for one electron
 summ = 0;
 for (i = 0; i < NEQ_vib; i++) {
  summ += pow(ydata[i],2);
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

 // Compute the analytical population on the c states
 Analytical_c(tout, numOutputSteps, energy, k_bandedge, k_bandtop, y);

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
 for (i = 1; i <= numsteps; ++i) {
  t = (tout*((double) i)/((double) numsteps));
  flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
  cout << endl << "CVode flag at step " << i << ": " << flag << endl;
#endif
  if (i % (numsteps/numOutputSteps) == 0) {
   fprintf(stdout, "\r%-.2lf percent done", ((double)i/((double)numsteps))*100);
   Output_checkpoint(
#ifdef DEBUG
     realImaginary, 
#endif
     allprob, yout, t, tkprob, tlprob, tcprob, tbprob, vibprob, times, qd_est,
     qd_est_diag, energy_expectation, (i*numOutputSteps/numsteps), energy,
     k_bandedge, k_bandtop, k_pops);
  }
 }
#ifdef DEBUG
 fclose(realImaginary);
#endif

 // compute final outputs //
 Compute_final_outputs(allprob, times, tkprob,
   tlprob, tcprob, tbprob, vibprob, energy,
   energy_expectation, numOutputSteps, qd_est, qd_est_diag);
 
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
