#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <time.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>

//#define DEBUG				// DANGER! Only turn on DEBUGf for small test runs, 
//#define DEBUGf				// otherwise output is enormous
//#define DEBUG_SAI			// debuggery related to checking against Sai's code
using namespace std;

// GLOBAL VARIABLES GO HERE //
 void * cvode_mem;			// pointer to block of CVode memory
 realtype * user_data;
 N_Vector y, yout;			// arrays of populations
 int Nk;				// number of each type of state
 int Nc;
 int Nb;
 int N_vib;				// number of vibronic states
 int Ik;				// index starters for each type of state
 int Ic;
 int Ib;
 int Ik_vib;				// index starters for each type of vibronic state
 int Ic_vib;
 int Ib_vib;
 int NEQ;				// total number of states/equations
 int NEQ_vib;
 double E_vib;					// vibrational energy
 double gkc;					// g factor between k and c states
 double gkb;					// g factor between k and b states
 double gbc;					// g factor between b and c states
 double gbb;					// g factor between b states
 realtype ** V;				// pointer to k-c coupling constants
 realtype ** FCkc;			// Franck-Condon factors
 realtype ** FCkb;
 realtype ** FCbb;
 realtype ** FCbc;
 realtype * energy;
 realtype * Vbridge;			// pointer to array of bridge coupling constants.
 					// first element [0] is Vkb1, last [Nb] is VcbN
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


void Build_k_energies(realtype * kEnergies, int numberOfKStates, realtype kBandEdge, realtype kBandTop) {
 
 int i;

 kEnergies[0] = kBandEdge;	// the bottom of the conduction band is set
 
 // loop over the remaining states.  This way the top of the band will be at kBandTop
 for (i = 1; i < numberOfKStates; i++) {
  kEnergies[i] = kEnergies[i-1] + (kBandTop-kBandEdge)/(numberOfKStates-1);
 }
}


void Build_k_pops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp) {

 int i;

 for (i = 0; i < Nk; i++)
  kPops[i] = 1.0/(1.0 + exp((kEnergies[i]-kBandEdge+0.1)*3.185e5/(temp)));	// dunno where the actual Fermi level is
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
 int scaleV = 1;
 realtype Vkc = 0.007349968763;

 if ((scaleV == 1) && (Nk > 1))
  Vkc = Vkc/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);

 for (i = 0; i < dim; i++)			// initialize
  for (j = 0; j < dim; j++)
   vArray[i][j] = 0.0;
 if (Nb == 0)					// no bridge
  // Vkc
  for (i = 0; i < Nk; i++)
   for (j = 0; j < Nc; j++)
    vArray[Ik+i][Ic+j] = Vkc;
 if (Nb > 0) {					// bridge
  // coupling between k and b1
  if ((scaleV == 1) && (Nk > 1)) {
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = Vbridge[0]/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);
    vArray[Ib][Ik+i] = Vbridge[0]/sqrt(Nk-1)*sqrt((kBandTop-kBandEdge)*27.211);
   }
  }
  else {
   for (i = 0; i < Nk; i++) {
    vArray[Ik+i][Ib] = Vbridge[0];
    vArray[Ib][Ik+i] = Vbridge[0];
   }
  }
  // coupling between bN and c
  if ((scaleV == 1) && (Nc > 1)) {
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
#ifdef DEBUG
 cout << "\nCoupling matrix:\n";
 for (i = 0; i < dim; i++) {
  for (j = 0; j < dim; j++)
   cout << vArray[i][j] << " ";
  cout << endl;
 }
#endif
 
}


// gives f(y,t) //
int f(realtype t, N_Vector y, N_Vector ydot, void * data) {
 
 int i, j, n, m;
 int IkRe, IkIm, IcRe, IcIm;
 realtype sinn;
 realtype coss;
 realtype Vee;
 // initialization
 for (i = 0; i < 2*NEQ_vib; i++)
  NV_Ith_S(ydot, i) = 0;

 if (Nb == 0) {						// no bridge
  realtype Vkc = V[Ik][Ic];
  for (i = 0; i < Nk; i++)
   for (j = 0; j < Nc; j++)
    for (n = 0; n < N_vib; n++)
     for (m = 0; m < N_vib; m++) {
      IkRe = Ik_vib + i*N_vib + n;			// indices
      IkIm = IkRe + NEQ_vib;
      IcRe = Ic_vib + j*N_vib + m;
      IcIm = IcRe + NEQ_vib;
      coss = cos((energy[IkRe] - energy[IcRe])*t);	// do the cosine now for speed
      sinn = sin((energy[IkRe] - energy[IcRe])*t);	// do the sine now for speed
      /* Franck-Condon indices should be consistent in direction from one end of the 
       * system to the other.  That is, if the first index is for the beginning state
       * the second should be for the next in the chain, and so on. */
      Vee = Vkc*FCkc[n][m];
      NV_Ith_S(ydot, IkRe) += Vee*(coss*NV_Ith_S(y, IcIm) + sinn*NV_Ith_S(y, IcRe)); // k Re
      NV_Ith_S(ydot, IkIm) += Vee*(sinn*NV_Ith_S(y, IcIm) - coss*NV_Ith_S(y, IcRe)); // k Im
      NV_Ith_S(ydot, IcRe) += Vee*(coss*NV_Ith_S(y, IkIm) - sinn*NV_Ith_S(y, IkRe)); // k Re
      NV_Ith_S(ydot, IcIm) -= Vee*(sinn*NV_Ith_S(y, IkIm) + coss*NV_Ith_S(y, IkRe)); // k Im
#ifdef DEBUGf
      cout << endl << "IkRe " << IkRe << " IkIm " << IkIm << " IcRe " << IcRe << " IcIm ";
      cout << IcIm << " V " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
#endif
     }
 }
 if (Nb > 0) {						// bridge
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
     coss = cos((energy[IkRe] - energy[IbRe])*t);
     sinn = sin((energy[IkRe] - energy[IbRe])*t);
     Vee = Vkb*FCkb[n][m];
     NV_Ith_S(ydot, IkRe) += Vee*(coss*NV_Ith_S(y, IbIm) + sinn*NV_Ith_S(y, IbRe));	// k Re
     NV_Ith_S(ydot, IkIm) += Vee*(sinn*NV_Ith_S(y, IbIm) - coss*NV_Ith_S(y, IbRe));	// k Im
     NV_Ith_S(ydot, IbRe) += Vee*(coss*NV_Ith_S(y, IkIm) - sinn*NV_Ith_S(y, IkRe));	// b1 Re
     NV_Ith_S(ydot, IbIm) -= Vee*(sinn*NV_Ith_S(y, IkIm) + coss*NV_Ith_S(y, IkRe));	// b1 Im
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
     coss = cos((energy[IcRe] - energy[IbRe])*t);
     sinn = sin((energy[IcRe] - energy[IbRe])*t);
     Vee = Vbc*FCbc[m][n];
     NV_Ith_S(ydot, IcRe) += Vee*(coss*NV_Ith_S(y, IbIm) + sinn*NV_Ith_S(y, IbRe));	// c Re
     NV_Ith_S(ydot, IcIm) += Vee*(sinn*NV_Ith_S(y, IbIm) - coss*NV_Ith_S(y, IbRe));	// c Im
     NV_Ith_S(ydot, IbRe) += Vee*(coss*NV_Ith_S(y, IcIm) - sinn*NV_Ith_S(y, IcRe));	// b1 Re
     NV_Ith_S(ydot, IbIm) -= Vee*(sinn*NV_Ith_S(y, IcIm) + coss*NV_Ith_S(y, IcRe));	// b1 Im
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
     coss = cos((energy[IbRe] - energy[IBRe])*t);
     sinn = sin((energy[IbRe] - energy[IBRe])*t);
     Vee = Vbridge[i+1]*FCbb[n][m];
     NV_Ith_S(ydot, IbRe) += Vee*(coss*NV_Ith_S(y, IBIm) + sinn*NV_Ith_S(y, IBRe));	// b Re
     NV_Ith_S(ydot, IbIm) += Vee*(sinn*NV_Ith_S(y, IBIm) - coss*NV_Ith_S(y, IBRe));	// b Im
     NV_Ith_S(ydot, IBRe) += Vee*(coss*NV_Ith_S(y, IbIm) - sinn*NV_Ith_S(y, IbRe));	// b1 Re
     NV_Ith_S(ydot, IBIm) -= Vee*(sinn*NV_Ith_S(y, IbIm) + coss*NV_Ith_S(y, IbRe));	// b1 Im
#ifdef DEBUGf
     cout << endl << "IbRe " << IbRe << " IbIm " << IbIm << " IBRe " << IBRe << " IBIm ";
     cout << IBIm << " Vbb " << Vee << " cos " << coss << " sin " << sinn << " t " << t;
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
    fprintf(dydt, "Re(c%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ic_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(c%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ic_vib + i*N_vib + j + NEQ_vib)*41.3414);
   }
  for (i = 0; i < Nb ; i++)
   for (j = 0; j < N_vib ; j++) {
    fprintf(dydt, "Re(b%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ib_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(b%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ib_vib + i*N_vib + j + NEQ_vib)*41.3414);
   }
  for (i = 0; i < Nk ; i++)
   for (j = 0; j < N_vib ; j++) {
    fprintf(dydt, "Re(k%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ik_vib + i*N_vib + j)*41.3414);
    fprintf(dydt, "Im(k%d,%-.5lf): %-.9lf\n", j, t, NV_Ith_S(ydot, Ik_vib + i*N_vib + j + NEQ_vib)*41.3414);
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


int Output_checkpoint(FILE * outputFile, FILE * kprobFile, N_Vector outputData, realtype time,
  realtype * totK, realtype * totC, realtype ** totB, realtype ** vibProb, realtype * times,
  realtype * energy_expectation, int index, realtype * energies) {

 int i, j, k, l;				// counters
 int Idx;
 realtype sumkpop = 0;
 realtype sumcpop = 0;
 realtype temp;

 for (i = 0; i < Nk; i++) {			// k populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ik_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  fprintf(outputFile, "%-.7lf ", temp);
  fprintf(kprobFile, "%-.7lf %-.7lf %-.7lf\n", time, energies[Ik_vib + i*N_vib], temp);
  sumkpop += temp;
 }
 fprintf(kprobFile, "\n");			// makes a blank line for gnuplot
 for (i = 0; i < Nc; i++) {			// c populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ic_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  fprintf(outputFile, "%-.7lf ", temp);
  sumcpop += temp;
 }
 for (i = 0; i < Nb; i++) {			// b populations
  temp = 0.0;
  for (j = 0; j < N_vib; j++) {
   Idx = Ib_vib + i*N_vib + j;			// indices
   temp += pow(NV_Ith_S(outputData, Idx),2) + pow(NV_Ith_S(outputData, Idx+NEQ_vib),2);
  }
  fprintf(outputFile, "%-.7lf ", temp);
  totB[i][index] = temp;
 }
 fprintf(outputFile, "\n");
 for (i = 0; i < N_vib; i++) {
  temp = 0;
  for (j = 0; j < NEQ; j++) {
   temp += pow(NV_Ith_S(outputData,i+j*N_vib),2) + pow(NV_Ith_S(outputData, i+j*N_vib+NEQ_vib),2);
  }
  vibProb[index][i] = temp;
 }
 totK[index] = sumkpop;
 totC[index] = sumcpop;
 times[index] = time;

 // temp += energies[i]*(pow(NV_Ith_S(outputData,i),2) + pow(NV_Ith_S(outputData,i+NEQ_vib),2));

 temp = 0.0;
 for (i = 0; i < NEQ; i++) {		// loop over all states
  for (j = 0; j < N_vib; j++) {		// loop over all vibronic states
   for (k = 0; k < NEQ; k++) {		// loop over all states (for coupling)
    for (l = 0; l < N_vib; l++) {	// loop over all vibronic states (for coupling)
    }
   }
  }
 }
 energy_expectation[index] = temp;

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


void Compute_final_outputs (realtype * time, realtype * tk, realtype * tc, realtype ** tb,
 realtype ** vibProb, realtype * energies, realtype * energy_expectation, int num) {

 FILE * tkprob;
 FILE * tcprob;
 FILE * vibprob;
 FILE * Ikprob;
 FILE * Icprob;
 FILE * kmax;
 FILE * cmax;
 FILE * totprob;
 FILE * energy;
 FILE * times;
 FILE * energy_exp;
 int i, j;
 realtype summ;

 tkprob = fopen("tkprob.out", "w");
 tcprob = fopen("tcprob.out", "w");
 vibprob = fopen("vibprob.out", "w");
 totprob = fopen("totprob.out", "w");
 Ikprob = fopen("Ikprob.out", "w");
 Icprob = fopen("Icprob.out", "w");
 kmax = fopen("kmax.out", "w");
 cmax = fopen("cmax.out", "w");
 energy = fopen("energy.out", "w");
 times = fopen("times.out", "w");
 energy_exp = fopen("energy_exp.out", "w");

 if (Nb > 0) {
  FILE * tbprob;
  FILE * Ibprob;
  FILE * bmax;

  tbprob = fopen("tbprob.out", "w");
  Ibprob = fopen("Ibprob.out", "w");
  bmax = fopen("bmax.out", "w");

  for (i = 0; i <= num; i++) {
   fprintf(tbprob, "%-.7lf %-.7lf", time[i], tb[0][i]);
   for (j = 1; j < Nb; j++)
    fprintf(tbprob, " %-.7lf", tb[j][i]);
   fprintf(tbprob, "\n");
  }
  
  fprintf(Ibprob, "%-.7lf", Integrate_arrays(tb[0], time, num+1));
  for (i = 1; i < Nb; i++)
   fprintf(Ibprob, " %-.7lf", Integrate_arrays(tb[i], time, num+1));
  
  fprintf(bmax, "%-.7lf", Find_array_maximum(tb[0], num+1));
  for (i = 1; i < Nb; i++)
  
   fprintf(bmax, "%-.7lf", Find_array_maximum(tb[i], num+1));
  fclose(tbprob);
  fclose(Ibprob);
  fclose(bmax);
 }

 fprintf(Ikprob, "%-.7lf", Integrate_arrays(tk, time, num+1));
 
 for (i = 0; i <= num; i++) {
  fprintf(tkprob, "%-.7lf %-.7lf\n", time[i], tk[i]);
  fprintf(tcprob, "%-.7lf %-.7lf\n", time[i], tc[i]);
  summ = tk[i] + tc[i];
  if (Nb > 0)
   for (j = 0; j < Nb; j++)
    summ += tb[j][i];
  fprintf(totprob, "%-.7lf %-.15lf\n", time[i], summ);
 }

 for (i = 0; i <= num; i++) {
  fprintf(vibprob, "%-.7lf %-.7lf", time[i], vibProb[i][0]);
  for (j = 1; j < N_vib; j++)
   fprintf(vibprob, " %-.7lf", vibProb[i][j]);
  fprintf(vibprob, "\n");
 }

 fprintf(Icprob, "%-.7lf", Integrate_arrays(tc, time, num+1));

 fprintf(kmax, "%-.7lf", Find_array_maximum(tk, num+1));
 fprintf(cmax, "%-.7lf", Find_array_maximum(tc, num+1));

 for (i = 0; i < Nk; i++)
  fprintf(energy, "%-.7lf\n", energies[Ik_vib + i*N_vib]);
 for (i = 0; i < Nc; i++)
  fprintf(energy, "%-.7lf\n", energies[Ic_vib + i*N_vib]);
 for (i = 0; i < Nb; i++) {
  fprintf(energy, "%-.7lf ", energies[Ib_vib + i*N_vib]);
  if (i < (Nb - 1))
   fprintf(energy, "\n");
 }

 for (i = 0; i <= num; i++)
  fprintf(times, "%-.7lf\n", time[i]);

 for (i = 0; i <= num; i++)
  fprintf(energy_exp, "%-.7lf\n", energy_expectation[i]);

 fclose(tkprob);
 fclose(tcprob);
 fclose(vibprob);
 fclose(totprob);
 fclose(Ikprob);
 fclose(Icprob);
 fclose(kmax);
 fclose(cmax);
 fclose(energy);
 fclose(times);
}

int main (int argc, char * argv[]) {

 // VARIABLES GO HERE//
 int i, j;					// counter!
 int flag;
 realtype * k_pops;				// pointers to arrays of populations
 realtype * c_pops;
 realtype * b_pops;
 realtype * ydata;				// pointer to ydata (contains all populations)
 realtype * k_energies;				// pointers to arrays of energies
 realtype * c_energies;
 realtype * b_energies;
 realtype k_bandedge;				// lower edge of bulk conduction band
 realtype k_bandtop;				// upper edge of bulk conduction band
 realtype temperature;				// system temperature
 realtype t0 = 0.0;				// initial time
 realtype t = 0;
 realtype tret;					// time returned by the solver
 time_t startRun;				// time at start of log
 time_t endRun;					// time at end of log
 struct tm * currentTime;			// time structure for localtime
 FILE * output;
 FILE * kprob;					// file formatted for plotting in gnuplot
 FILE * log;
 realtype * tkprob; 				// total probability in k, c, b states at each timestep
 realtype * tcprob;
 realtype ** tbprob;
 realtype ** vibprob;
 realtype * times;
 realtype * energy_expectation;			// expectation value of energy at each timestep
 // END VARIABLES //

 // OPEN LOG FILE; PUT IN START TIME //
 log = fopen("log.out", "w");			// note that this file is closed at the end of the program
 time(&startRun);
 currentTime = localtime(&startRun);
 fprintf(log, "Run started at %s\n", asctime(currentTime));
 
 // ASSIGN VARIABLES FROM RUN SCRIPT //
 realtype abstol = atof(argv[1]);		// absolute tolerance (for SUNDIALS)
 realtype reltol = atof(argv[2]);		// relative tolerance (for SUNDIALS)
 realtype tout = atof(argv[3]);			// final time reached by solver in atomic units
 int numsteps = atoi(argv[4]);			// number of time steps
 int numOutputSteps = atoi(argv[5]);
 // bulk parameters //
 k_bandedge = atof(argv[6]);			// lower band edge of conduction band
 k_bandtop = atof(argv[7]);			// upper band edge of bulk conduction band
 Nk = atoi(argv[8]);				// number of k states
 // physical parameters //
 temperature = atof(argv[9]);			// temperature of the system
 // vibronic parameters //
 N_vib = atoi(argv[10]);			// number of vibronic states
 E_vib = atof(argv[11]);			// vibrational energy
 gkc = atof(argv[12]);				// g factor between k and c states
 gkb = atof(argv[13]);				// g factor between k and b states
 gbc = atof(argv[14]);				// g factor between b and c states
 gbb = atof(argv[15]);				// g factor between b states
#ifdef DEBUG
 cout << endl;
 for (i = 0; i < argc; i++)
  cout << "argv[" << i << "] is " << argv[i] << endl;
#endif
 // DONE ASSIGNING VARS FROM SCRIPT //


 // READ DATA FROM INPUTS //
 Nc = Number_of_values("ins/c_energies.in");
 Nb = Number_of_values("ins/b_energies.in");
 k_pops = new realtype [Nk];
 c_pops = new realtype [Nc];
 b_pops = new realtype [Nb];
 k_energies = new realtype [Nk];
 c_energies = new realtype [Nc];
 b_energies = new realtype [Nb];
 Vbridge = new realtype [Nb+1];
 Read_array_from_file(c_pops, "ins/c_pops.in", Nc);
 Read_array_from_file(c_energies, "ins/c_energies.in", Nc);
 if ( Nb > 0) {
  Read_array_from_file(b_energies, "ins/b_energies.in", Nb);
  Read_array_from_file(Vbridge, "ins/Vbridge.in", Nb + 1);
 }
 // DONE READING //
#ifdef DEBUG
 cout << "\nDone reading things from inputs.\n";
#endif

 // PREPROCESS DATA FROM INPUTS //
 tkprob = new realtype [numOutputSteps+1];	// total population on k, b, c at each timestep
 tcprob = new realtype [numOutputSteps+1];
 tbprob = new realtype * [numOutputSteps+1];
 for (i = 0; i < Nb; i++)
  tbprob[i] = new realtype [numOutputSteps+1];
 vibprob = new realtype * [numOutputSteps+1];
 for (i = 0; i < numOutputSteps+1; i++)
  vibprob[i] = new realtype [N_vib];
 times = new realtype [numOutputSteps+1];
 energy_expectation = new realtype [numOutputSteps+1];	// expectation value of energy; for sanity checking
 NEQ = Nk+Nc+Nb;				// total number of equations set
 NEQ_vib = (Nk+Nc+Nb)*N_vib;
 Ik = 0;					// set index start positions for each type of state
 Ic = Nk;
 Ib = Ic+Nc;
 Ik_vib = 0;
 Ic_vib = Nk*N_vib;
 Ib_vib = Ic_vib + Nc*N_vib;
 Build_k_energies(k_energies, Nk, k_bandedge, k_bandtop);	// create bulk conduction quasicontinuum
 Initialize_array(b_pops, Nb, 0.0);		// populate b states
 Initialize_array(k_pops, Nk, 0.0);		// populate k states (all zero to start off)
 //Build_k_pops(k_pops, k_energies, k_bandedge, temperature);	// populate k states (all zero to start off)
 V = new realtype * [NEQ];
 for (i = 0; i < NEQ; i++)
  V[i] = new realtype [NEQ];
 Build_v(V, NEQ, k_bandedge, k_bandtop);	// assign k-c coupling constants
 if (Nb == 0) {					// assign Franck-Condon factors
  FCkc = new realtype * [N_vib];
  for (i = 0; i < N_vib; i++)
   FCkc[i] = new realtype [N_vib];
  Build_Franck_Condon_factors(FCkc, gkc, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCkc:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7lf ", FCkc[i][j]);
    cout << endl;
   }
#endif
 }
 else if (Nb > 0) {
  if (Nb > 1) {
   FCbb = new realtype * [N_vib];
   for (i = 0; i < N_vib; i++)
    FCbb[i] = new realtype [N_vib];
   Build_Franck_Condon_factors(FCbb, gbb, N_vib, N_vib);
#ifdef DEBUG
   cout << "\n FCbb:\n";
   for (i = 0; i < N_vib; i++) {
    for (j = 0; j < N_vib; j++)
     printf("%+.7lf ", FCbb[i][j]);
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
     printf("%+.7lf ", FCkb[i][j]);
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
     printf("%+.7lf ", FCbc[i][j]);
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
 for (i = 0; i < Nk; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[k(" << i << "," << j << ")] = " << ydata[Ik_vib + i*N_vib + j + NEQ_vib] << endl;
 for (i = 0; i < Nc; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[c(" << i << "," << j << ")] = " << ydata[Ic_vib + i*N_vib + j + NEQ_vib] << endl;
 for (i = 0; i < Nb; i++)
  for (j = 0; j < N_vib; j++)
  cout << "starting ydata: Im[b(" << i << "," << j << ")] = " << ydata[Ib_vib + i*N_vib + j + NEQ_vib] << endl;
 cout << endl;
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
 cout << endl;
#endif
 // DONE PREPROCESSING //

 // Creates N_Vector y with initial populations which will be used by CVode//
 y = N_VMake_Serial(2*NEQ_vib, ydata);
 yout = N_VClone(y);

 // print t = 0 information //
 Normalize_NV(y, 1.00);			// normalizes all populations to 1; this is for one electron
 output = fopen("allprob.out", "w");	// note that output is closed after advancing the solution in time
 kprob = fopen("kprob.out", "w");
 Output_checkpoint(output, kprob, y, t0, tkprob, tcprob, tbprob, vibprob, times, energy_expectation, 0, energy);

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
  if (i % (numsteps/numOutputSteps) == 0)
   Output_checkpoint(output, kprob, yout, t, tkprob, tcprob, tbprob, vibprob, times, energy_expectation, (i*numOutputSteps/numsteps), energy);
 }
 fclose(output);				// note that file is opened when printing t=0 information
 fclose(kprob);

 // compute final outputs //
 Compute_final_outputs(times, tkprob, tcprob, tbprob, vibprob, energy, energy_expectation, numOutputSteps);
 
 // finalize log file //
 time(&endRun);
 currentTime = localtime(&endRun);
 fprintf(log, "Run ended at %s\n", asctime(currentTime));
 fprintf(log, "Run took %.3lf seconds.", difftime(endRun, startRun));
 fclose(log);					// note that the log file is opened after variable declaration
 printf("\nRun took %.3lf seconds.\n", difftime(endRun, startRun));

 // deallocate memory for N_Vectors //
 N_VDestroy_Serial(y);
 N_VDestroy_Serial(yout);

 // free solver memory //
 CVodeFree(&cvode_mem);

 // delete all these guys
 delete [] tkprob;
 delete [] tcprob;
 delete [] tbprob;
 delete [] vibprob;
 delete [] k_pops;
 delete [] c_pops;
 delete [] b_pops;
 delete [] energy;
 delete [] V;
 delete [] Vbridge;
 delete [] k_energies;
 delete [] c_energies;
 delete [] b_energies;
 delete [] ydata;
 delete [] FCkc;
 delete [] FCkb;
 delete [] FCbc;
 delete [] FCbb;
 fprintf(stderr, "\nwhoo\n");
 return 0;
}
