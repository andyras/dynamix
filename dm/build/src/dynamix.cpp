#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <time.h>
#include <numeric>
#include <complex>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>
#include <map>
#include <fftw/fftw3.h>

#include "libdynamix_input_parser.h"
#include "libdynamix_outputs.h"
#include "output.h"
#include "numerical.h"

/* DEBUG compiler flag: turn on to generate basic debug outputs.         */
#define DEBUG
// DEBUG2 flag: turn on for more numerical output
//#define DEBUG2
/* DANGER! Only turn on DEBUGf for small test runs, otherwise output is       */
/* enormous (many GB).  This flag turns on debug output within the f          */
/* function.                                                                  */
//#define DEBUGf

using namespace std;

// GLOBAL VARIABLES GO HERE //
void * cvode_mem;			// pointer to block of CVode memory
realtype * user_data;
N_Vector y, yout;			// arrays of populations
int Nk;				// number of each type of state
int Nc;
int Nb;
int Nl;
int Ik;				// index starters for each type of state
int Ic;
int Ib;
int Il;
int NEQ;				// total number of states/equations
int NEQ2;
int numOutputSteps;			// number of timesteps
realtype k_bandedge;			// lower edge of bulk conduction band
realtype k_bandtop;			// upper edge of bulk conduction band
double muLK;                           // transition dipole moment from l to k (energy a.u.)
double pumpFWHM;                       // FWHM of pump pulse (time a.u.)
double pumpPeak;                       // time of peak of pump pulse (a.u.)
double pumpFreq;                       // frequency of pump pulse (energy a.u.)
double pumpAmpl;                       // intensity of pump pulse (electric field a.u.)
double pumpPhase;                      // pump pulse phase (in units of radians)
realtype ** V;				// pointer to k-c coupling constants
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

void buildCoupling (realtype ** vArray, int dim, realtype kBandEdge,
                    realtype kBandTop, realtype * energy,
		    std::map<std::string, bool> &outs) {
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

 if (outs["couplings.out"]) {
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
 
}


int f(realtype t, N_Vector y, N_Vector ydot, void * data) {
// gives f(y,t) for CVODE

 // TODO: decompose 'data' variable into energies
 // TODO: only calculate upper/lower triangle of diagonal blocks in DM
 
 //// Checklist for debugging this code
 // * Make sure all indices have declared variables
 // * Check that indices are computed properly
 //   * I_{x_i}Re = I_x + i;
 //   * I_{x_i,y_j}Re = NEQ*I_{x_i} + I_{y_j};
 //   * I_{x_i,y_j}Im = I_{x_i,y_j}Re + NEQ2;
 // * Check that indices used have correct Re/Im
 // * Check signs of each term and each sum
 //   * In general Re/Im have opposite signs
 //   * For complex conjugates, Re have same sign, Im opposite
 // * Check usage of y (LHS) and ydot (RHS)
 // * Check usage of sumRe vs. sumIm
 
 // variables for indices
 int IkRe, IkIm, IbRe, IbIm, IcRe, IcIm;
 int IkpRe, IbpRe, IcpRe;
 int IijRe, IijIm, IjiRe, IjiIm;
 int IkkpRe, IkkpIm, IbbpRe, IbbpIm, IccpRe, IccpIm;
 int IcpcRe, IcpcIm;
 int IkbRe, IkbIm, IkcRe, IkcIm;
 int IbkRe, IbkIm, IbcRe, IbcIm;
 int IbkpRe, IbkpIm, IbcpRe, IbcpIm;
 int IckRe, IckIm, IcbRe, IcbIm;

 // accumulators
 realtype sumRe, sumIm;

 // energy gap
 realtype wij;

 // initialize ydot
 for (int ii = 0; ii < 2*NEQ2; ii++) {
  NV_Ith_S(ydot, ii) = 0;
 }

 // equations of motion for system with bridge
 if (bridge_on) {
  // couplings
  realtype Vkb = V[Ik][Ib];
  realtype Vbc = V[Ib][Ic];
  // These indices won't change with only one bridge
  IbRe = Ib;
  IbbpRe = NEQ*IbRe + IbRe;
  IbbpIm = IbbpRe + NEQ2;

  //// Contributions from energy gaps
  for (int ii = 0; ii < NEQ; ii++) {
   // only need off-diagonal terms (would be 0 if ii == jj)
   for (int jj = 0; jj < ii; jj++) {
    IijRe = NEQ*ii + jj;
    IijIm = IijRe + NEQ2;
    IjiRe = NEQ*jj + ii;
    IjiIm = IjiRe + NEQ2;
    wij = energy[ii] - energy[jj];

    // real parts of complex conjugates are the same
    NV_Ith_S(ydot, IijRe) += wij*NV_Ith_S(y, IijIm);
    NV_Ith_S(ydot, IjiRe) += wij*NV_Ith_S(y, IijIm);
    // imaginary parts of complex conjugates have opposite sign
    NV_Ith_S(ydot, IijIm) += wij*NV_Ith_S(y, IijRe);
    NV_Ith_S(ydot, IjiIm) -= wij*NV_Ith_S(y, IijRe);
   }
  }

  // \dot{\rho_{kk'}}
  sumRe = 0.0;
  sumIm = 0.0;
  for (int ii = 0; ii < Nk; ii++) {
   IkRe = Ik + ii;
   for (int jj = 0; jj < Nk; jj++) {
    IkpRe = Ik + jj;
    IkkpRe = NEQ*IkRe + IkpRe;
    IkkpIm = IkkpRe + NEQ2;
    IbkpRe = NEQ*IbRe + IkpRe;
    IbkpIm = IbkpRe + NEQ2;
    IkbRe = NEQ*Ik + IbRe;
    IkbIm = IkbRe + NEQ2;

    // real part	V_{kb}(\rho_{bk'} - \rho_{kb}) term
    sumRe += NV_Ith_S(y, IbkpIm) - NV_Ith_S(y, IkbIm);
    // imaginary part	V_{kb}(\rho_{bk'} - \rho_{kb}) term
    sumIm -= NV_Ith_S(y, IbkpRe) - NV_Ith_S(y, IkbRe);
   }
  }

  // real part		V_{kb}(\rho_{bk'} - \rho_{kb}) term
  NV_Ith_S(ydot, IkkpRe) += Vkb*sumRe;
  // imaginary part	V_{kb}(\rho_{bk'} - \rho_{kb}) term
  NV_Ith_S(ydot, IkkpIm) += Vkb*sumIm;

  //// \dot{\rho_{bb}}
  sumRe = 0.0;
  sumIm = 0.0;
  /// contribution from \rho_{bk}, \rho_{kb}
  for (int ii = 0; ii < Nk; ii++) {
   IkRe = Ik + ii;
   IkbRe = NEQ*IkRe + IbRe;
   IkbIm = IkbRe + NEQ2;
   IbkRe = NEQ*IbRe + IkRe;
   IbkIm = IbkRe + NEQ2;

   // real part		V_{kb}(\dot{\rho_{kb}} - \dot{\rho_{bk}}) term
   sumRe += NV_Ith_S(y, IkbIm) - NV_Ith_S(y, IbkIm);
   // imaginary part	V_{kb}(\dot{\rho_{kb}} - \dot{\rho_{bk}}) term
   sumIm -= NV_Ith_S(y, IkbRe) - NV_Ith_S(y, IbkRe);
  }

  // real part		V_{kb}(\dot{\rho_{kb}} - \dot{\rho_{bk}}) term
  NV_Ith_S(ydot, IbbpRe) += Vkb*sumRe;
  // imaginary part	V_{kb}(\dot{\rho_{kb}} - \dot{\rho_{bk}}) term
  NV_Ith_S(ydot, IbbpIm) += Vkb*sumIm;

  /// contribution from \rho_{bc}, \rho_{cb}
  for (int ii = 0; ii < Nc; ii++) {
   IcRe = Ic + ii;
   IcbRe = NEQ*IcRe + IbRe;
   IcbIm = IcbRe + NEQ2;
   IbcRe = NEQ*IbRe + IcRe;
   IbcIm = IbcRe + NEQ2;

   // real part		V_{bc}(\dot{\rho_{bc}} - \dot{\rho_{cb}}) term
   sumRe += NV_Ith_S(y, IbcIm) - NV_Ith_S(y, IcbIm);
   // imaginary part	V_{bc}(\dot{\rho_{bc}} - \dot{\rho_{cb}}) term
   sumIm -= NV_Ith_S(y, IbcRe) - NV_Ith_S(y, IcbRe);
  }

  // real part		V_{bc}(\dot{\rho_{bc}} - \dot{\rho_{cb}}) term
  NV_Ith_S(ydot, IbbpRe) += Vbc*sumRe;
  // imaginary part	V_{bc}(\dot{\rho_{bc}} - \dot{\rho_{cb}}) term
  NV_Ith_S(ydot, IbbpIm) += Vbc*sumIm;

  //// \dot{\rho_{cc'}}
  sumRe = 0.0;
  sumIm = 0.0;
  for (int ii = 0; ii < Nc; ii++) {
   IcRe = Ic + ii;
   for (int jj = 0; jj < Nc; jj++) {
    IcpRe = Ic + jj;
    IccpRe = NEQ*IcRe + IcpRe;
    IccpIm = IccpRe + NEQ2;
    IbcpRe = NEQ*IbRe + IcpRe;
    IbcpIm = IbcpRe + NEQ2;
    IcbRe = NEQ*IcRe + IbRe;
    IcbIm = IcbRe + NEQ2;

    // real part	V_{bc}(\rho_{bc'} - \rho_{cb}) term
    sumRe += NV_Ith_S(y, IbcpIm) - NV_Ith_S(y, IcbIm);
    // imaginary part	V_{bc}(\rho_{bc'} - \rho_{cb}) term
    sumIm -= NV_Ith_S(y, IbcpRe) - NV_Ith_S(y, IcbRe);
   }
  }

  // real part		V_{bc}(\rho_{bc'} - \rho_{cb}) term
  NV_Ith_S(ydot, IccpRe) += Vbc*sumRe;
  // imaginary part	V_{bc}(\rho_{bc'} - \rho_{cb}) term
  NV_Ith_S(ydot, IccpIm) += Vbc*sumIm;

  //// \dot{\rho_{kb}} and \dot{\rho_{bk}}
  for (int ii = 0; ii < Nk; ii++) {
   IkRe = Ik + ii;
   IkbRe = NEQ*IkRe + IbRe;
   IkbIm = IkbRe + NEQ2;
   IbkRe = NEQ*IbRe + IkRe;
   IbkIm = IbkRe + NEQ2;

   // real part		V_{kb}\rho_{bb} term
   sumRe = NV_Ith_S(y, IbbpIm);
   // imaginary part	V_{kb}\rho_{bb} term
   sumIm = -1*NV_Ith_S(y, IbbpRe);

   for (int jj = 0; jj < Nk; jj++) {
    IkpRe = Ik + jj;
    IkkpRe = NEQ*IkRe + IkpRe;
    IkkpIm = IkkpRe + NEQ2;
    // real part	V_{kb}\rho_{kk'} term
    sumRe -= NV_Ith_S(y, IkkpIm);
    // imaginary part	V_{kb}\rho_{kk'} term
    sumIm += NV_Ith_S(y, IkkpRe);
   }
   // real part		V_{kb} terms
   NV_Ith_S(ydot, IkbRe) += Vkb*sumRe;
   NV_Ith_S(ydot, IbkRe) += Vkb*sumRe;
   // imaginary part	V_{kb} terms
   NV_Ith_S(ydot, IkbIm) += Vkb*sumIm;
   NV_Ith_S(ydot, IbkIm) -= Vkb*sumIm;

   sumRe = 0.0;
   sumIm = 0.0;
   for (int jj = 0; jj < Nc; jj++) {
    IcRe = Ic + jj;
    IkcRe = NEQ*IkRe + IcRe;
    IkcIm = IkcRe + NEQ2;

    // real part	V_{bc}\rho_{kc} term
    sumRe -= NV_Ith_S(y, IkcIm);
    // imaginary part	V_{bc}\rho_{kc} term
    sumIm += NV_Ith_S(y, IkcRe);
   }

   // real part		V_{bc}\rho_{kc} term
   NV_Ith_S(ydot, IkbRe) += Vbc*sumRe;
   NV_Ith_S(ydot, IbkRe) += Vbc*sumRe;
   // imaginary part	V_{bc}\rho_{kc} term
   NV_Ith_S(ydot, IkbIm) += Vbc*sumIm;
   NV_Ith_S(ydot, IbkIm) += Vbc*sumIm;
  }

  //// \dot{\rho_{kc}} and \dot{\rho_{ck}}
  for (int ii = 0; ii < Nk; ii++) {
   IkRe = Ik + ii;
   for (int jj = 0; jj < Nc; jj++) {
    IcRe = Ic + jj;
    IkcRe = NEQ*IkRe + IcRe;
    IkcIm = IkcRe + NEQ2;
    IbcRe = NEQ*IbRe + IcRe;
    IbcIm = IbcRe + NEQ2;
    IkbRe = NEQ*IkRe + IbRe;
    IkbIm = IkbRe + NEQ2;

    // real part	V_{kb}\rho_{bc} term
    NV_Ith_S(ydot, IkcRe) += Vkb*NV_Ith_S(y, IbcIm) - Vbc*NV_Ith_S(y, IkbIm);
    NV_Ith_S(ydot, IckRe) += Vkb*NV_Ith_S(y, IbcIm) - Vbc*NV_Ith_S(y, IkbIm);
    // imaginary part	V_{kb}\rho_{bc} term
    NV_Ith_S(ydot, IkcIm) -= Vkb*NV_Ith_S(y, IbcRe) - Vbc*NV_Ith_S(y, IkbRe);
    NV_Ith_S(ydot, IckIm) += Vkb*NV_Ith_S(y, IbcRe) - Vbc*NV_Ith_S(y, IkbRe);
   }
  }

  //// \dot{\rho_{bc}} and \dot{\rho_{cb}}
  for (int ii = 0; ii < Nc; ii++) {
   IcRe = Ic + ii;
   IcIm = IcRe + NEQ;
   IbcRe = NEQ*IbRe + IcRe;
   IbcIm = IbcRe + NEQ2;
   IcbRe = NEQ*IcRe + IbRe;
   IcbIm = IcbRe + NEQ2;

   // real part		V_bc\rho_{bb} term
   sumRe = -1*NV_Ith_S(y, IbbpIm);
   // imaginary part	V_bc\rho_{bb} term
   sumIm = NV_Ith_S(y, IbbpRe);

   for (int jj = 0; jj < Nc; jj++) {
    IcpRe = Ic + jj;
    IcpcRe = NEQ*IcpRe + IcRe;
    IcpcIm = IcpcRe + NEQ2;

    // real part	V_bc\rho_{c'c} term
    sumRe += NV_Ith_S(y, IcpcIm);
    // imaginary part	V_bc\rho_{c'c} term
    sumIm -= NV_Ith_S(y, IcpcRe);
   }

   // real part		V_bc\rho_{c'c} term
   NV_Ith_S(ydot, IbcRe) += Vbc*sumRe;
   NV_Ith_S(ydot, IcbRe) += Vbc*sumRe;
   // imaginary part	V_bc\rho_{c'c} term
   NV_Ith_S(ydot, IbcIm) += Vbc*sumIm;
   NV_Ith_S(ydot, IcbIm) -= Vbc*sumIm;

   sumRe = 0.0;
   sumIm = 0.0;
   for (int jj = 0; jj < Nk; jj++) {
    IkRe = Ik + jj;
    IkcRe = NEQ*IkRe + IcRe;
    IkcIm = IkcRe + NEQ2;
    // real part	V_kb\rho_{kc} term
    sumRe += NV_Ith_S(y, IkcIm);
    // imaginary part	V_kb\rho_{kc} term
    sumIm -= NV_Ith_S(y, IkcRe);
   }

    // real part	V_kb\rho_{kc} term
    NV_Ith_S(ydot, IbcRe) += Vkb*sumRe;
    NV_Ith_S(ydot, IcbRe) += Vkb*sumRe;
    // imaginary part	V_kb\rho_{kc} term
    NV_Ith_S(ydot, IbcIm) += Vkb*sumIm;
    NV_Ith_S(ydot, IcbIm) -= Vkb*sumIm;
  }
 }

 return 0;
}


/*
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
 for (i = 0; i < NEQ2; i++) {
  fprintf(realImaginary, "%-.9e %-.9e\n", NV_Ith_S(outputData, i), NV_Ith_S(outputData, i+NEQ2));
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
*/

void Compute_final_outputs (double ** allprobs, realtype * time, realtype * tk,
  realtype * tl, realtype * tc, realtype * tb, realtype * energies,
  realtype * energy_expectation, int num, double * qd_est, double * qd_est_diag,
  std::map<std::string, bool> &outs) {
// makes output files for time-dependent properties

 FILE * tkprob;
 FILE * tlprob;
 FILE * tcprob;
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

 if (outs["tkprob.out"]) {
  tkprob = fopen("tkprob.out", "w");
 }
 if (outs["tlprob.out"]) {
  tlprob = fopen("tlprob.out", "w");
 }
 if (outs["tcprob.out"]) {
  tcprob = fopen("tcprob.out", "w");
 }
 if (outs["kprobs.out"]) {
  kprobs = fopen("kprobs.out", "w");
 }
 if (outs["kprobs_gnuplot.out"]) {
  kprobs_gnuplot = fopen("kprobs_gnuplot.out", "w");
 }
 if (outs["cprobs_gnuplot.out"]) {
  cprobs_gnuplot = fopen("cprobs_gnuplot.out", "w");
 }
 if (outs["cprobs.out"]) {
  cprobs = fopen("cprobs.out", "w");
 }
 if (outs["bprobs.out"]) {
  bprobs = fopen("bprobs.out", "w");
 }
 if (outs["totprob.out"]) {
  totprob = fopen("totprob.out", "w");
 }
 if (outs["Ikprob.out"]) {
  Ikprob = fopen("Ikprob.out", "w");
 }
 if (outs["Icprob.out"]) {
  Icprob = fopen("Icprob.out", "w");
 }
 if (outs["kmax.out"]) {
  kmax = fopen("kmax.out", "w");
 }
 if (outs["cmax.out"]) {
  cmax = fopen("cmax.out", "w");
 }
 if (outs["cmax_t.out"]) {
  cmax_t = fopen("cmax_t.out", "w");
 }
 if (outs["cmax_first.out"]) {
  cmax_first = fopen("cmax_first.out", "w");
 }
 if (outs["cmax_first_t.out"]) {
  cmax_first_t = fopen("cmax_first_t.out", "w");
 }
 if (outs["energy.out"]) {
  energy = fopen("energy.out", "w");
 }
 if (outs["times.out"]) {
  times = fopen("times.out", "w");
 }
 if (outs["pump_intensity.out"]) {
  pump_intensity = fopen("pump_intensity.out", "w");
 }
 if (outs["energy_exp.out"]) {
  energy_exp = fopen("energy_exp.out", "w");
 }
 if (outs["qd_est.out"]) {
  qd_estimate = fopen("qd_est.out", "w");
 }
 if (outs["qd_est_diag.out"]) {
  qd_estimate_diag = fopen("qd_est_diag.out", "w");
 }

 if (outs["kprobs.out"]) {
  for (i = 0 ; i < num ; i++) {				// print k probabilities over time
   fprintf(kprobs, "%-.7g", time[i]);
   for (j = 0; j < Nk ; j++)
    fprintf(kprobs, " %-.7g", allprobs[i][Ik+j]);
   fprintf(kprobs, "\n");
  }
 }

 if (outs["kprobs_gnuplot.out"]) {
  for (i = 0 ; i < num ; i++) {
   for (j = 0 ; j < Nk ; j++ )
    fprintf(kprobs_gnuplot, "%-.7g %-.7g %-.7g\n", time[i], energies[Ik + j], allprobs[i][Ik+j]);
   if (i < (num - 1))
    fprintf(kprobs_gnuplot, "\n");			// makes a blank line for gnuplot
  }
 }

 if (outs["cprobs.out"]) {
  for (i = 0 ; i < num ; i++) {				// print c probabilities over time
   fprintf(cprobs, "%-.7g", time[i]);
   for (j = 0; j < Nc ; j++)
    fprintf(cprobs, " %-.7g", allprobs[i][Ic+j]);
   fprintf(cprobs, "\n");
  }
 }

 if (outs["cprobs_gnuplot.out"]) {
  for (i = 0 ; i < num ; i++) {
   for (j = 0 ; j < Nc ; j++ )
    fprintf(cprobs_gnuplot, "%-.7g %-.7g %-.7g\n", time[i], energies[Ic + j], allprobs[i][Ic+j]);
   if (i < (num - 1))
    fprintf(cprobs_gnuplot, "\n");			// makes a blank line for gnuplot
  }
 }

 if (outs["bprobs.out"]) {
  for (i = 0 ; i < num ; i++) {				// print b probabilities over time
   fprintf(bprobs, "%-.7g", time[i]);
   for (j = 0; j < Nb ; j++)
    fprintf(bprobs, " %-.7g", allprobs[i][Ib+j]);
   fprintf(bprobs, "\n");
  }
 }

 if (Nb > 0) {
  FILE * tbprob;
  FILE * Ibprob;
  FILE * bmax;
  // TODO find out why this is set to 0 and then written
  double max_b_prob = 0;

  if (outs["tbprob.out"]) {
   tbprob = fopen("tbprob.out", "w");
   for (i = 0; i <= num; i++)				// print total b population
    fprintf(tbprob, "%-.7g %-.7g\n", time[i], tb[i]);
   fclose(tbprob);
  }

  if (outs["Ibprob.out"]) {
   Ibprob = fopen("Ibprob.out", "w");
   fprintf(Ibprob, "%-.7g", integrateArray(tb, time, num+1));
   fclose(Ibprob);
  }

  if (outs["bmax.out"]) {
   bmax = fopen("bmax.out", "w");
   for (i = 0 ; i < num + 1 ; i++) {
    for (j = 0 ; j < Nb ; j++) {
     if (allprobs[i][Ib+j] > max_b_prob) {
      max_b_prob = allprobs[i][Ib+j];
     }
    }
   }
   // TODO THE MYSTERY!!!
   //fprintf(bmax, "%-.7g", findArrayMaximum(tb, num+1));
   fprintf(bmax, "%-.7g", max_b_prob);
   fclose(bmax);
  }
 }

 for (i = 0; i <= num; i++) {				// print total k, c population
  if (outs["tkprob.out"]) {
   fprintf(tkprob, "%-.7g %-.7g\n", time[i], tk[i]);
  }
  if (outs["tlprob.out"]) {
   fprintf(tlprob, "%-.7g %-.7g\n", time[i], tl[i]);
  }
  if (outs["tcprob.out"]) {
   fprintf(tcprob, "%-.7g %-.7g\n", time[i], tc[i]);
  }
  summ = tk[i] + tc[i] + tl[i];
  if (Nb > 0)
   summ += tb[i];
  if (outs["totprob.out"]) {
   fprintf(totprob, "%-.7g %-.15g\n", time[i], summ);
  }
 }

 if (outs["Ikprob.out"]) {
  fprintf(Ikprob, "%-.7g", integrateArray(tk, time, num+1));
 }
 if (outs["Icprob.out"]) {
  fprintf(Icprob, "%-.7g", integrateArray(tc, time, num+1));
 }
 
 if (outs["kmax.out"]) {
  fprintf(kmax, "%-.7g", findArrayMaximum(tk, num+1));
 }
 if (outs["cmax.out"]) {
  fprintf(cmax, "%-.7g", findArrayMaximum(tc, num+1));
 }
 if (outs["cmax_first.out"]) {
  fprintf(cmax_first, "%-.7g", findFirstArrayMaximum(tc, num+1));
 }
 if (outs["cmax_t.out"]) {
  fprintf(cmax_t, "%-.7g", time[findArrayMaximumIndex(tc, num+1)]);
 }
 if (outs["cmax_first_t.out"]) {
  fprintf(cmax_first_t, "%-.7g", time[findFirstArrayMaximumIndex(tc, num+1)]);
 }

 if (outs["energy.out"]) {
  // energy.out should be all the energies on one row, since it's used for
  // the movie maker.
  fprintf(energy, "%-.7g", energies[0]);
  for (i = 1; i < NEQ; i++)
   fprintf(energy, " %-.7g", energies[i]);
 }

 for (i = 0; i <= num; i++) {
  if (outs["times.out"]) {
   fprintf(times, "%-.7g\n", time[i]);
  }
  if (outs["pump_intensity.out"]) {
   fprintf(pump_intensity, "%-.7g %-.7g\n", time[i], pump(time[i], pumpFWHM, pumpAmpl, pumpPeak, pumpFreq, pumpPhase));
  }
  if (outs["energy_exp.out"]) {
   fprintf(energy_exp, "%-.7g %-.7g\n", time[i], energy_expectation[i]);
  }
  if (outs["qd_estimate.out"]) {
   fprintf(qd_estimate, "%-.7g %-.7g\n", time[i], qd_est[i]);
  }
  if (outs["qd_estimate_diag.out"]) {
   fprintf(qd_estimate_diag, "%-.7g %-.7g\n", time[i], qd_est_diag[i]);
  }
 }

 if (outs["tkprob.out"]) {
  fclose(tkprob);
 }
 if (outs["tlprob.out"]) {
  fclose(tlprob);
 }
 if (outs["tcprob.out"]) {
  fclose(tcprob);
 }
 if (outs["kprobs.out"]) {
  fclose(kprobs);
 }
 if (outs["kprobs_gnuplot.out"]) {
  fclose(kprobs_gnuplot);
 }
 if (outs["cprobs.out"]) {
  fclose(cprobs);
 }
 if (outs["cprobs_gnuplot.out"]) {
  fclose(cprobs_gnuplot);
 }
 if (outs["bprobs.out"]) {
  fclose(bprobs);
 }
 if (outs["totprob.out"]) {
  fclose(totprob);
 }
 if (outs["Ikprob.out"]) {
  fclose(Ikprob);
 }
 if (outs["Icprob.out"]) {
  fclose(Icprob);
 }
 if (outs["kmax.out"]) {
  fclose(kmax);
 }
 if (outs["cmax.out"]) {
  fclose(cmax);
 }
 if (outs["cmax_first.out"]) {
  fclose(cmax_first);
 }
 if (outs["cmax_t.out"]) {
  fclose(cmax_t);
 }
 if (outs["cmax_first_t.out"]) {
  fclose(cmax_first_t);
 }
 if (outs["energy.out"]) {
  fclose(energy);
 }
 if (outs["times.out"]) {
  fclose(times);
 }
 if (outs["qd_estimate.out"]) {
  fclose(qd_estimate);
 }
 if (outs["qd_estimate_diag.out"]) {
  fclose(qd_estimate_diag);
 }

 // compute derivatives
 FILE * tkDeriv;
 double * tkderivs = new double [numOutputSteps-5];
 Derivative(tk, numOutputSteps, tkderivs, time[1]-time[0]);
 if (outs["tkderiv.out"]) {
  tkDeriv = fopen("tkderiv.out","w");
  for (i = 0; i < numOutputSteps-5; i++) {
   fprintf(tkDeriv, "%-.7g %-.9g\n", time[i+2], tkderivs[i]);
  }
  fclose(tkDeriv);
 }
 delete [] tkderivs;

 FILE * tcDeriv;
 double * tcderivs = new double [numOutputSteps-5];
 Derivative(tc, numOutputSteps, tcderivs, time[1]-time[0]);
 if (outs["tcderiv.out"]) {
  tcDeriv = fopen("tcderiv.out","w");
  for (i = 0; i < numOutputSteps-5; i++) {
   fprintf(tcDeriv, "%-.7g %-.9g\n", time[i+2], tcderivs[i]);
  }
  fclose(tcDeriv);
 }
 delete [] tcderivs;

 // compute rates.  Rate is calculated as (dP(t)/dt)/P(t), where P(t) is
 // the population at time t.
 FILE * tkRate;
 if (outs["tkrate.out"]) {
  tkRate = fopen("tkrate.out","w");
  for (i = 0; i < numOutputSteps-5; i++) {
   fprintf(tkRate, "%-.7g %-.9g\n", time[i+2], tkderivs[i]/tk[i+2]);
  }
  fclose(tkRate);
 }

 FILE * tcRate;
 if (outs["tcrate.out"]) {
  tcRate = fopen("tcrate.out","w");
  for (i = 0; i < numOutputSteps-5; i++) {
   fprintf(tcRate, "%-.7g %-.9g\n", time[i+2], tcderivs[i]/tc[i+2]);
  }
  fclose(tcRate);
 }

}


/* Makes a gnuplot file to plot the QD populations over time */
void plot_cprobs(int n, double t, double k_bandtop, double k_bandedge, int Nk) {
 std::ofstream output("cprobs.plt");
 output << "#!/usr/bin/env gnuplot\n\n"
 << "reset\n"
 << "set terminal pdfcairo enhanced size 4in,3in font 'Arial-Bold,14'\n"
 << "set output '/dev/null'\n"
 << "!transpose -o _transpose ../outs/cprobs.out\n"
 << "plot '../outs/cprobs_transpose.out' every :::1 u ($1*" << t << "/" << n << "):(-$2):3 matrix with image\n"
 << "set output 'cprobs.pdf'\n"
 << "set title 'Electron probability density in QD'\n"
 << "set border 0\n"
 << "unset ytics\n"
 << "set xtics scale 0\n"
 << "set ylabel 'States above band edge'\n"
 << "set xlabel 'Time (a.u.)'\n"
 << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]\n"
 << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]\n"
 << "unset key\n"
 << "unset colorbox\n"
 << "set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee', 4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n"
 << "repl\n";

 return;
}

/* Updates \rho(t) at each time step. */
void updateDM(N_Vector dm, realtype * dmt, int timeStep) {
 for (int ii = 0; ii < NEQ2; ii++) {
  dmt[NEQ2*timeStep + ii] = NV_Ith_S(dm, ii);
  dmt[NEQ2*timeStep + ii + NEQ2] = NV_Ith_S(dm, ii + NEQ2);
 }

 return;
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
 realtype * wavefunction;			// (initial) wavefunction
 realtype * dm;					// density matrix
 realtype * dmt;				// density matrix in time
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
 double ** allprob;				// populations in all states at all times
 realtype * times;
 realtype * qd_est;
 realtype * qd_est_diag;
 realtype * energy_expectation;			// expectation value of energy at each timestep
 const char * inputFile;			// name of input file
 inputFile = "ins/parameters.in";
 std::map<std::string, bool> outs;		// map of output file names to bool
 // END VARIABLES //
 
 // Decide which output files to make
#ifdef DEBUG
 std::cout << "Assigning outputs as specified in " << inputFile << "\n";
#endif
 assignOutputs(inputFile, outs);
#ifdef DEBUG
 // print out which outputs will be made
 for (std::map<std::string, bool>::iterator it = outs.begin(); it != outs.end(); it++) {
  std::cout << "Output file: " << it->first << " will be created.\n";
 }
#endif

 // OPEN LOG FILE; PUT IN START TIME //
 if (outs["log.out"]) {
  log = fopen("log.out", "w");			// note that this file is closed at the end of the program
 }
 time(&startRun);
 currentTime = localtime(&startRun);
 if (outs["log.out"]) {
  fprintf(log, "Run started at %s\n", asctime(currentTime));
 }
 
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
 // laser parameters
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

 if (outs["log.out"]) {
  // make a note about the laser intensity.
  fprintf(log,"The laser intensity is %.5e W/cm^2.\n\n",pow(pumpAmpl,2)*3.5094452e16);
 }

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
 Nc = numberOfValuesInFile("ins/c_energies.in");
 Nb = numberOfValuesInFile("ins/b_energies.in");
 k_pops = new realtype [Nk];
 c_pops = new realtype [Nc];
 b_pops = new realtype [Nb];
 l_pops = new realtype [Nl];
 k_energies = new realtype [Nk];
 c_energies = new realtype [Nc];
 b_energies = new realtype [Nb];
 l_energies = new realtype [Nl];
 if (numberOfValuesInFile("ins/c_pops.in") != Nc) {
  fprintf(stderr, "ERROR [Inputs]: c_pops and c_energies not the same length.\n");
  return -1;
 }
 readArrayFromFile(c_energies, "ins/c_energies.in", Nc);
 if (bridge_on) {
  if (bridge_on && (Nb < 1)) {
   cerr << "\nERROR: bridge_on but no bridge states.  The file b_energies.in is probably empty.\n";
   return -1;
  }
  Vbridge = new realtype [Nb+1];
  readArrayFromFile(b_energies, "ins/b_energies.in", Nb);
  readArrayFromFile(Vbridge, "ins/Vbridge.in", Nb + 1);
 }
 else {
  Nb = 0;
  Vnobridge = new realtype [1];
  readArrayFromFile(Vnobridge, "ins/Vnobridge.in", 1);
 }
 // DONE READING //
#ifdef DEBUG
 cout << "\nDone reading things from inputs.\n";
#endif

 // PREPROCESS DATA FROM INPUTS //
 NEQ = Nk+Nc+Nb+Nl;				// total number of equations set
 NEQ2 = NEQ*NEQ;				// number of elements in DM
#ifdef DEBUG
 cout << "\nTotal number of states: " << NEQ << endl;
 cout << Nk << " bulk, " << Nc << " QD, " << Nb << " bridge, " << Nl << " bulk VB.\n";
#endif
 tkprob = new realtype [numOutputSteps+1];	// total population on k, b, c at each timestep
 tcprob = new realtype [numOutputSteps+1];
 tbprob = new realtype [numOutputSteps+1];
 tlprob = new realtype [numOutputSteps+1];
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
 // assign bulk conduction band energies
 buildContinuum(k_energies, Nk, k_bandedge, k_bandtop);
 // assign bulk valence band energies
 buildContinuum(l_energies, Nl, k_bandedge - valenceBand - bulk_gap, k_bandedge - bulk_gap);

 //// Build initial wavefunction

 // bridge states (empty to start)
 initializeArray(b_pops, Nb, 0.0);

 // coefficients in bulk and other states depend on input conditions in bulk
 if (bulk_constant) {
#ifdef DEBUG
  cout << "\ninitializing k_pops\n";
#endif
  initializeArray(k_pops, Nk, 0.0);
#ifdef DEBUG
  cout << "\ninitializing k_pops with constant probability in range of states\n";
#endif
  initializeArray(k_pops+Nk_first-1, Nk_final-Nk_first+1, 1.0);
#ifdef DEBUG
  cout << "\nThis is k_pops:\n";
  for (i = 0; i < Nk; i++) {
   cout << k_pops[i] << endl;
  }
  cout << "\n";
#endif
  initializeArray(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  initializeArray(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (bulk_Gauss) {
  buildKPopsGaussian(k_pops, k_energies, k_bandedge,
                        bulkGaussSigma, bulkGaussMu, Nk);   // populate k states with FDD
  initializeArray(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  initializeArray(c_pops, Nc, 0.0);		// QD states empty to start
 }
 else if (qd_pops) {
  readArrayFromFile(c_pops, "ins/c_pops.in", Nc);	// QD populations from file
  initializeArray(l_pops, Nl, 0.0);		// populate l states (all 0 to start off)
  initializeArray(k_pops, Nk, 0.0);             // populate k states (all zero to start off)
 }
 else {
  initializeArray(k_pops, Nk, 0.0);             // populate k states (all zero to start off)
  initializeArray(l_pops, Nl, 1.0);		// populate l states (all populated to start off)
  initializeArray(c_pops, Nc, 0.0);		// QD states empty to start
 }

 // create empty wavefunction
 wavefunction = new realtype [2*NEQ];
 initializeArray(wavefunction, 2*NEQ, 0.0);

 // assign real parts of wavefunction coefficients (imaginary are zero)
 for (i = 0; i < Nk; i++)
  wavefunction[Ik + i] = k_pops[i];
 for (i = 0; i < Nc; i++)
  wavefunction[Ic + i] = c_pops[i];
 for (i = 0; i < Nb; i++)
  wavefunction[Ib + i] = b_pops[i];
 for (i = 0; i < Nl; i++)
  wavefunction[Il + i] = l_pops[i];

 if (outs["psi_start.out"]) {
  outputWavefunction(wavefunction, NEQ);
 }

 // Give all coefficients a random phase
 if (random_phase) {
  float phi;
  // set the seed
  if (random_seed == -1) { srand(time(NULL)); }
  else { srand(random_seed); }
  for (i = 0; i < NEQ; i++) {
   phi = 2*3.1415926535*(float)rand()/(float)RAND_MAX;
   wavefunction[i] = wavefunction[i]*cos(phi);
   wavefunction[i + NEQ] = wavefunction[i + NEQ]*sin(phi);
  }
 }

#ifdef DEBUG
 // print out details of wavefunction coefficients
 cout << endl;
 for (i = 0; i < Nk; i++)
  cout << "starting wavefunction: Re[k(" << i << ")] = " << wavefunction[Ik + i] << endl;
 for (i = 0; i < Nc; i++)
  cout << "starting wavefunction: Re[c(" << i << ")] = " << wavefunction[Ic + i] << endl;
 for (i = 0; i < Nb; i++)
  cout << "starting wavefunction: Re[b(" << i << ")] = " << wavefunction[Ib + i] << endl;
 for (i = 0; i < Nl; i++)
  cout << "starting wavefunction: Re[l(" << i << ")] = " << wavefunction[Il + i] << endl;
 for (i = 0; i < Nk; i++)
  cout << "starting wavefunction: Im[k(" << i << ")] = " << wavefunction[Ik + i + NEQ] << endl;
 for (i = 0; i < Nc; i++)
  cout << "starting wavefunction: Im[c(" << i << ")] = " << wavefunction[Ic + i + NEQ] << endl;
 for (i = 0; i < Nb; i++)
  cout << "starting wavefunction: Im[b(" << i << ")] = " << wavefunction[Ib + i + NEQ] << endl;
 for (i = 0; i < Nl; i++)
  cout << "starting wavefunction: Im[l(" << i << ")] = " << wavefunction[Il + i + NEQ] << endl;
 cout << endl;
 summ = 0;
 for (i = 0; i < 2*NEQ; i++) {
  summ += pow(wavefunction[i],2);
 }
 cout << "\nTotal population is " << summ << "\n\n";
#endif

 // Assemble array of energies
 energy = new realtype [NEQ];
 for (i = 0; i < Nk; i++)
  energy[Ik + i] = k_energies[i];
 for (i = 0; i < Nc; i++)
  energy[Ic + i] = c_energies[i];
 for (i = 0; i < Nb; i++)
  energy[Ib + i] = b_energies[i];
 for (i = 0; i < Nl; i++)
  energy[Il + i] = l_energies[i];
 user_data = energy;

#ifdef DEBUG
 // output energies
 for (i = 0; i < Nk; i++)
  cout << "energy[k(" << i << ")] = " << energy[Ik + i] << endl;
 for (i = 0; i < Nc; i++)
  cout << "energy[c(" << i << ")] = " << energy[Ic + i] << endl;
 for (i = 0; i < Nb; i++)
  cout << "energy[b(" << i << ")] = " << energy[Ib + i] << endl;
 for (i = 0; i < Nl; i++)
  cout << "energy[l(" << i << ")] = " << energy[Il + i] << endl;
 cout << endl;
#endif

 // assign coupling constants
 V = new realtype * [NEQ];
 for (i = 0; i < NEQ; i++)
  V[i] = new realtype [NEQ];
 buildCoupling(V, NEQ, k_bandedge, k_bandtop, energy, outs);

 if (outs["log.out"]) {
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
 }

 // Create the initial density matrix
 dm = new realtype [2*NEQ2];
 initializeArray(dm, 2*NEQ2, 0.0);
 for (int ii = 0; ii < NEQ; ii++) {
  // diagonal part
  dm[NEQ*ii + ii] = pow(wavefunction[ii],2) + pow(wavefunction[ii + NEQ],2);
  // off-diagonal part
  for (int jj = 0; jj < ii; jj++) {
   // real part of \rho_{ii,jj}
   dm[NEQ*ii + jj] = wavefunction[ii]*wavefunction[jj] + wavefunction[ii+NEQ]*wavefunction[jj+NEQ];
   // imaginary part of \rho_{ii,jj}
   dm[NEQ*ii + jj + NEQ2] = wavefunction[ii]*wavefunction[jj+NEQ] - wavefunction[jj]*wavefunction[ii*NEQ];
   // real part of \rho_{jj,ii}
   dm[NEQ*jj + ii] = dm[NEQ*ii + jj];
   // imaginary part of \rho_{jj,ii}
   dm[NEQ*jj + ii + NEQ2] = -1*dm[NEQ*ii + jj + NEQ*NEQ];
  }
 }

 // Create the array to store the density matrix in time
 dmt = new realtype [2*NEQ2*numOutputSteps];
 initializeArray(dmt, 2*NEQ2*numOutputSteps, 0.0);

#ifdef DEBUG
 // print out density matrix
 cout << "\nDensity matrix before normalization:\n\n";
 for (int ii = 0; ii < NEQ; ii++) {
  for (int jj = 0; jj < NEQ; jj++) {
   fprintf(stdout, "(%+.1e,%+.1e) ", dm[NEQ*ii + jj], dm[NEQ*ii + jj + NEQ2]);
  }
  fprintf(stdout, "\n");
 }
#endif

 // Normalize the DM so that populations add up to 1.
 summ = 0.0;
 for (int ii = 0; ii < NEQ; ii++) {
  // assume here that diagonal elements are all real
  summ += dm[NEQ*ii + ii];
 }
 if ( summ == 0.0 ) {
  cerr << "\nFATAL ERROR [populations]: total population is 0!\n";
  return -1;
 }
 if (summ != 1.0) {
  // the variable 'summ' is now a multiplicative normalization factor
  summ = 1.0/summ;
  for (int ii = 0; ii < 2*NEQ2; ii++) {
   dm[ii] *= summ;
  }
 }
#ifdef DEBUG
 cout << "\nThe normalization factor for the density matrix is " << summ << "\n\n";
#endif

 // Error checking for total population; recount population first
 summ = 0.0;
 for (int ii = 0; ii < NEQ; ii++) {
  summ += dm[NEQ*ii + ii];
 }
 if ( fabs(summ-1.0) > 1e-12 ) {
  cerr << "\nWARNING [populations]: After normalization, total population is not 1, it is " << summ << "!\n";
 }
#ifdef DEBUG
 cout << "\nAfter normalization, the sum of the populations in the density matrix is " << summ << "\n\n";
#endif

 // DONE PREPROCESSING //

 // Creates N_Vector y with initial populations which will be used by CVode//
 y = N_VMake_Serial(2*NEQ2, dm);
 // put in t = 0 information
 updateDM(y, dmt, 0);
 // the vector yout has the same dimensions as y
 yout = N_VClone(y);

 // print t = 0 information //
 summ = 0;
 for (i = 0; i < 2*NEQ2; i++) {
  summ += pow(NV_Ith_S(y, i),2);
 }
#ifdef DEBUG
 realImaginary = fopen("real_imaginary.out", "w");
#endif
 /*
 Output_checkpoint(
#ifdef DEBUG
   realImaginary, 
#endif
   allprob, y, t0, tkprob, tlprob, tcprob, tbprob, times, qd_est,
   qd_est_diag, energy_expectation, 0, energy, k_bandedge, k_bandtop, k_pops);
   */

 // create CVode object
 // this is a stiff problem, I guess?
 cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
 // make 'energy' array available to CVode via 'user_data'
 flag = CVodeSetUserData(cvode_mem, (void *) user_data);

 // initialize CVode solver //
 flag = CVodeInit(cvode_mem, &f, t0, y);

 // specify integration tolerances //
 flag = CVodeSStolerances(cvode_mem, reltol, abstol);

 // attach linear solver module //
 flag = CVDense(cvode_mem, 2*NEQ2);

 // advance the solution in time! //
 // use CVODE for time-dependent H
 for (i = 1; i <= numsteps; ++i) {
  t = (tout*((double) i)/((double) numsteps));
  flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
  cout << endl << "CVode flag at step " << i << ": " << flag << endl;
#endif
  if (i % (numsteps/numOutputSteps) == 0) {
   fprintf(stderr, "\r%-.2lf percent done", ((double)i/((double)numsteps))*100);
   updateDM(yout, dmt, i*numOutputSteps/numsteps);
   /*
   Output_checkpoint(
#ifdef DEBUG
     realImaginary, 
#endif
     allprob, yout, t, tkprob, tlprob, tcprob, tbprob, times, qd_est,
     qd_est_diag, energy_expectation, (i*numOutputSteps/numsteps), energy,
     k_bandedge, k_bandtop, k_pops);
   */
  }
 }

 // compute final outputs //
 Compute_final_outputs(allprob, times, tkprob,
   tlprob, tcprob, tbprob, energy,
   energy_expectation, numOutputSteps, qd_est, qd_est_diag, outs);
 computeDMOutput(dmt, NEQ, V, energy, times, numOutputSteps, outs);

 outputDMt(dmt, NEQ, numOutputSteps, outs);

 // compute time-independent outputs
 if (outs["energy.out"]) {
  FILE * energyFile = fopen("energy.out", "w");
  for (i = 0; i < NEQ; i++) {
   fprintf(energyFile, "%-.9e\n", energy[i]);
  }
  fclose(energyFile);
 }
#ifdef DEBUG
 fclose(realImaginary);
#endif

 if (outs["cprobs.plt"] && (Nc > 1)) {
  plot_cprobs(numOutputSteps, tout, k_bandtop, k_bandedge, Nk);
 }
 
 // finalize log file //
 time(&endRun);
 currentTime = localtime(&endRun);
 if (outs["log.out"]) {
  fprintf(log, "Final status of 'flag' variable: %d\n\n", flag);
  fprintf(log, "Run ended at %s\n", asctime(currentTime));
  fprintf(log, "Run took %.3g seconds.\n", difftime(endRun, startRun));
  fclose(log);					// note that the log file is opened after variable declaration
 }
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
 fprintf(stderr, "\nwhoo\n");
 return 0;
}
