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
#include <omp.h>

#include "libdynamix_input_parser.h"
#include "libdynamix_outputs.h"
#include "output.h"
#include "numerical.h"
#include "params.h"
#include "userdata.h"

/* DEBUG compiler flag: turn on to generate basic debug outputs.         */
//#define DEBUG
// DEBUG2 flag: turn on for more numerical output
//#define DEBUG2
/* DANGER! Only turn on DEBUGf for small test runs, otherwise output is       */
/* enormous (many GB).  This flag turns on debug output within the f          */
/* function.                                                                  */
// #define DEBUGf
//#define DEBUG_DMf

using namespace std;

// Struct of parameters
PARAMETERS params;

// Hamiltonian
realtype * Ham;

// GLOBAL VARIABLES GO HERE //
#ifdef DEBUG_DMf
FILE * dmf;				// file for density matrix coeff derivatives in time
#endif
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

int f(realtype t, N_Vector y, N_Vector ydot, void * user_data) {
// gives f(y,t) for CVODE

 // user_data is a pointer to the params struct
 PARAMETERS * p;
 p = (PARAMETERS *) user_data;

 // extract parameters from p
 std::vector<realtype> H = p->H;
 int N = p->NEQ;
 int N2 = p->NEQ2;

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

#ifdef DEBUG_DMf
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


int main (int argc, char * argv[]) {

 // VARIABLES GO HERE//
 // number of processors
 int nproc;

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
 nproc = 0;
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
  else if (input_param == "nproc") { nproc = atof(param_val.c_str()); }
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
 cout << "nproc is " << nproc << endl;
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
 if (Nl < 0) {
  fprintf(stderr, "ERROR [Inputs]: Nl less than 0.\n");
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
  // feed coupling array to params
  params.Vbridge.resize(Nb+1);
  for (int ii = 0; ii < Nb + 1; ii++) {
   params.Vbridge[ii] = Vbridge[ii];
  }
 }
 else {
  Nb = 0;
  Vnobridge = new realtype [1];
  readArrayFromFile(Vnobridge, "ins/Vnobridge.in", 1);
  // feed coupling array to params
  params.Vnobridge.resize(1);
  params.Vnobridge[0] = Vbridge[0];
 }
 // DONE READING //
#ifdef DEBUG
 cout << "\nDone reading things from inputs.\n";
#endif

 // PREPROCESS DATA FROM INPUTS //
 // set number of processors for OpenMP
 omp_set_num_threads(nproc);

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
 for (i = 0; i <= numOutputSteps; i++) {
  allprob[i] = new double [NEQ];
 }
 // assign times.
 times = new realtype [numOutputSteps+1];
 for (int ii = 0; ii <= numOutputSteps; ii++) {
  times[ii] = float(ii)/numOutputSteps*tout;
 }
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

 //// feed parameters to struct
 params.Nk = Nk;
 params.Nc = Nc;
 params.Nb = Nb;
 params.Nl = Nl;
 params.Ik = Ik;
 params.Ic = Ic;
 params.Ib = Ib;
 params.Il = Il;
 params.NEQ = NEQ;
 params.NEQ2 = NEQ2;
 params.numOutputSteps = numOutputSteps;
 params.bridge_on = bridge_on;
 params.tout = tout;
 params.kBandEdge = k_bandedge;
 params.kBandTop = k_bandtop;
 params.scale_bubr = scale_bubr;
 params.scale_brqd = scale_brqd;
 params.scale_buqd = scale_buqd;
 params.parabolicCoupling = parabolicCoupling;
 params.times.resize(numOutputSteps);
 for (int ii = 0; ii <= numOutputSteps; ii++) {
  params.times[ii] = times[ii];
 }

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

 // feed energies to params
 params.energies.resize(NEQ);
 for (int ii = 0; ii < NEQ; ii++) {
  params.energies[ii] = energy[ii];
#ifdef DEBUG
  std::cout << "params.energies[" << ii << "] is " << params.energies[ii] << "\n";
#endif
 }

 // assign coupling constants
 V = new realtype * [NEQ];
 for (i = 0; i < NEQ; i++)
  V[i] = new realtype [NEQ];
 buildCoupling(V, &params, outs);

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
   dm[NEQ*ii + jj + NEQ2] = wavefunction[ii]*wavefunction[jj+NEQ] - wavefunction[jj]*wavefunction[ii+NEQ];
   // real part of \rho_{jj,ii}
   dm[NEQ*jj + ii] = dm[NEQ*ii + jj];
   // imaginary part of \rho_{jj,ii}
   dm[NEQ*jj + ii + NEQ2] = -1*dm[NEQ*ii + jj + NEQ*NEQ];
  }
 }

 // Create the array to store the density matrix in time
 dmt = new realtype [2*NEQ2*(numOutputSteps+1)];
 initializeArray(dmt, 2*NEQ2*(numOutputSteps+1), 0.0);

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

  // build Hamiltonian
#ifdef DEBUG
  fprintf(stderr, "Building Hamiltonian.\n");
#endif
  realtype * H = new realtype [NEQ2];
  buildHamiltonian(H, energy, V, &params);
  if (outs["ham.out"]) {
   outputSquareMatrix(H, NEQ, "ham.out");
  }
  Ham = new realtype [NEQ2];
  for (int ii = 0; ii < NEQ2; ii++) {
   Ham[ii] = H[ii];
  }
  // add Hamiltonian to params
  params.H.resize(NEQ2);
  for (int ii = 0; ii < NEQ2; ii++) {
   params.H[ii] = H[ii];
  }

 // DONE PREPROCESSING //
 
#ifdef DEBUG
 cout << "\nCreating N_Vectors.\n";
 cout << "\nProblem size is " << 2*NEQ2 << " elements.\n";
#endif
 // Creates N_Vector y with initial populations which will be used by CVode//
 y = N_VMake_Serial(2*NEQ2, dm);
 // put in t = 0 information
 updateDM(y, dmt, 0, &params);
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

 // create CVode object
 // this is a stiff problem, I guess?
#ifdef DEBUG
 cout << "\nCreating cvode_mem object.\n";
#endif
 cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
 flag = CVodeSetUserData(cvode_mem, (void *) &params);

#ifdef DEBUG
 cout << "\nInitializing CVode solver.\n";
#endif
 // initialize CVode solver //
 flag = CVodeInit(cvode_mem, &f, t0, y);

#ifdef DEBUG
 cout << "\nSpecifying integration tolerances.\n";
#endif
 // specify integration tolerances //
 flag = CVodeSStolerances(cvode_mem, reltol, abstol);

#ifdef DEBUG
 cout << "\nAttaching linear solver module.\n";
#endif
 // attach linear solver module //
 flag = CVDense(cvode_mem, 2*NEQ2);

 // advance the solution in time! //
 // use CVODE for time-dependent H
#ifdef DEBUG
 cout << "\nAdvancing the solution in time.\n";
#endif
#ifdef DEBUG_DMf
 cout << "Creating output file for density matrix coefficient derivatives in time.\n";
 dmf = fopen("dmf.out", "w");
#endif
 for (i = 1; i <= numsteps; ++i) {
  t = (tout*((double) i)/((double) numsteps));
  flag = CVode(cvode_mem, t, yout, &tret, 1);
#ifdef DEBUGf
  //cout << endl << "CVode flag at step " << i << ": " << flag << endl;
#endif
  if (i % (numsteps/numOutputSteps) == 0) {
   fprintf(stderr, "\r%-.2lf percent done", ((double)i/((double)numsteps))*100);
   updateDM(yout, dmt, i*numOutputSteps/numsteps, &params);
  }
 }
#ifdef DEBUG_DMf
 cout << "Closing output file for density matrix coefficients in time.\n";
 fclose(dmf);
#endif

#ifdef DEBUG
 fclose(realImaginary);
#endif

 if (outs["cprobs.plt"] && (Nc > 1)) {
  plot_cprobs(params);
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

 // Compute density matrix outputs.
#ifdef DEBUG
 cout << "Computing outputs...";
#endif
 computeDMOutput(dmt, outs, &params);
#ifdef DEBUG
 std::cout << "done.";
#endif

#ifdef DEBUG
 fprintf(stdout, "Deallocating N_Vectors.\n");
#endif
 //  TODO why does this block break?
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
