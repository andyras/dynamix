#include "numerical.hpp"
//#define DEBUG

//#define DEBUG_NUMERICAL
//#define DEBUG_READVECTOR

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

  fclose(input);

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

  return;
}

/* reads in values from a file to a vector. */
void readVectorFromFile(std::vector<realtype> & v, const char * fileName, int n) {
  // resize vector according to number of lines in file
  v.resize(numberOfValuesInFile(fileName));
#ifdef DEBUG_READVECTOR
  std::cout << numberOfValuesInFile(fileName) << " values in " << fileName << std::endl;
#endif
  if (numberOfValuesInFile(fileName) != n) {
    std::cerr << "ERROR reading in vector from " << fileName << ": wrong number of values in file." << std::endl;
    exit(1);
  }

  // read in the file
  std::ifstream in(fileName);

  for (int ii = 0; ii < n; ii++) {
    in >> v[ii];
  }

  in.close();

  return;
}

/* Returns an array of length n with all values set to initializeValue. */
void initializeArray(realtype * array, int n, realtype initializeValue) {
#ifdef DEBUG
  std::cout << "initializeValue is " << initializeValue << std::endl;
#endif

  for (int ii = 0; ii < n; ii++) {
    array[ii] = initializeValue;
  }

  return;
}

/* builds energies for a quasicontinuum (evenly spaced) */
void buildContinuum(realtype * Energies, int numberOfStates, realtype BandEdge, realtype BandTop) {
  if (numberOfStates == 0) {
    return;
  }

  Energies[0] = BandEdge;	// the bottom of the conduction band is set

  // loop over the remaining states.  This way the top of the band will be at BandTop
  for (int ii = 1; ii < numberOfStates; ii++) {
    Energies[ii] = Energies[ii-1] + (BandTop-BandEdge)/(numberOfStates-1);
  }

  return;
}

/* populates a set of states according to a Fermi-Dirac distribution.
 * I'm not sure where the actual Fermi level is, so it defaults to 0.01 Eh
 * below the lowest-energy state in the set of states being populated.
 */
void buildKPops(realtype * kPops, realtype * kEnergies, realtype kBandEdge, realtype temp, int Nk) {

  for (int ii = 0; ii < Nk; ii++) {
    kPops[ii] = sqrt(1.0/(1.0 + exp((kEnergies[ii]-kBandEdge+0.01)*3.185e5/(temp))));
#ifdef DEBUG
    std::cout << "\nk population at state " << ii << " is: "
      << sqrt(1.0/(1.0 + exp((kEnergies[ii]-kBandEdge+0.01)*3.185e5/(temp))));
#endif
  }
#ifdef DEBUG
  std::cout << std::endl;
#endif
}

/* populates a set of states according to a Gaussian distribution. */
void buildKPopsGaussian(realtype * kPops, realtype * kEnergies, realtype kBandEdge, double sigma, double mu, int Nk) {
#ifdef DEBUG_NUMERICAL
  std::cout << "Conduction band edge is " << kBandEdge << std::endl;
  std::cout << "Gaussian mu is          " << mu << std::endl;
  std::cout << "Gaussian sigma is       " << sigma << std::endl;
  std::cout << "Number of k states is   " << Nk << std::endl;
#endif

  for (int ii = 0; ii < Nk; ii++) {
    // take the square root so that populations have proper width given by sigma
    kPops[ii] = sqrt((1/(sigma*sqrt(2*3.1415926535)))*exp(-pow((kEnergies[ii]-(kBandEdge+mu)),2)/(2*pow(sigma,2))));
#ifdef DEBUG_NUMERICAL
    std::cout << "\nk population at state " << ii << " is: " << kPops[ii];
#endif
  }
#ifdef DEBUG_NUMERICAL
  std::cout << std::endl;
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

/* gives the value of a Gaussian laser pulse (electric field) at time t */
realtype gaussPulse(realtype t, double pumpFWHM, double pumpAmpl, double pumpPeak, double pumpFreq, double pumpPhase) {
  double sigma = pumpFWHM/2.35482005;	// conversion between FWHM and Gaussian 2nd moment
  return pumpAmpl*exp((-pow(t-pumpPeak, 2))/(2*pow(sigma, 2)))*cos(pumpFreq*t + pumpPhase);
}

/* normalizes an N_Vector so that the populations of all states are
 * normalized to value 'total'
 */
int Normalize_NV(N_Vector nv, realtype total) {
  realtype summ = 0;

  for (int ii = 0; ii < NV_LENGTH_S(nv); ii++) {
    summ += (NV_Ith_S(nv, ii)*NV_Ith_S(nv, ii));
  }
  summ = sqrt(summ);
  for (int ii = 0; ii < NV_LENGTH_S(nv); ii++) {
    NV_Ith_S(nv, ii) = total*NV_Ith_S(nv,ii)/summ;
  }

  return 0;
}

/* gives the value of [a + b*sin^2(c*t + d)] at a certain t */
double sin2(double a, double b, double c, double d, double t) {
  return a + b*pow(sin(c*t + d),2);
}

/* compute the six-point derivative of an array.
 * Assumes array elements are evenly spaced.
 * Output array is six elements shorter than input.
 */
int Derivative(double *inputArray, int inputLength, double *outputArray, double timestep) {
  if (inputLength < 6 ) {
    fprintf(stderr, "ERROR [Derivative]: array has length less than 6 elements, cannot proceed");
    return -1;
  }

  for (int ii = 2; ii < inputLength-3; ii++) {
    outputArray[ii-2] = (2* inputArray[ii+3]
	-15*inputArray[ii+2]
	+60*inputArray[ii+1]
	-20*inputArray[ii]
	-30*inputArray[ii-1]
	+3 *inputArray[ii-2])/(60*timestep);
  }

  return 0;
}

/* Compute the six-point time derivative of a 2D array.
 * Assumes array elements are evenly spaced in time.
 * Output array is six elements shorter than the input, skipping first two and
 * last three points.  The output array must be preallocated.
 */
void arrayDeriv(double * in, int nt, int m, int dim, double * out, double dt) {
  // n is length of array (number of time points)
  // m is width of array (number of populations)
  // dim is p->NEQ
  if (nt < 6) {
    std::cerr << "Error [" << __FUNCTION__ << "]: array must be at least six elements." << std::endl;
    exit(-1);
  }

  for (int ii = 0; ii < m; ii++) {
    for (int jj = 2; jj < (nt-3); jj++) {
      out[ii*nt + (jj-2)] = (2*in[ii*nt + jj+3]
	  - 15*in[ii*nt + jj+2]
	  + 60*in[ii*nt + jj+1]
	  - 20*in[ii*nt + jj]
	  - 30*in[ii*nt + jj-1]
	  +  3*in[ii*nt + jj-2])
	/
	(60*dt);
    }
  }

  return;
}

/* Riemann sum of an array (values) at time points (time).
 * Does not assume equal spacing in time.
 */
realtype integrateArray(realtype * values, realtype * time, int num) {
  realtype riemann = 0;

  for (int ii = 0; ii < num-1; ii++)
    riemann += (values[ii+1] + values[ii])*(time[ii+1]-time[ii])/2;

  return riemann;
}

/* Returns maximum element in an array. */
realtype findArrayMaximum(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];

  for (int ii = 1; ii < num; ii++)
    if (inputArray[ii] > currentMax)
      currentMax = inputArray[ii];

  return currentMax;
}

/* Finds the first maximum in an array (the first point where the next
 * point is smaller in value).
 */
realtype findFirstArrayMaximum(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax)
      currentMax = inputArray[ii];
    if (inputArray[ii] < currentMax)
      break;
  }

  return currentMax;
}

/* This function returns the index of the first maximum in an array.
 * Warning: will return 0 as index of first maximum if the second array element
 * is less than the first.  This may not be what you want.
 */
int findFirstArrayMaximumIndex(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];
  int currentMax_index = 0;

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax) {
      currentMax = inputArray[ii];
      currentMax_index = ii;
    }
    if (inputArray[ii] < currentMax)
      break;
  }

  return currentMax_index;
}

/* returns index of first maximum in an array. */
int findArrayMaximumIndex(realtype * inputArray, int num) {
  realtype currentMax = inputArray[0];
  int currentMax_index = 0;

  for (int ii = 1; ii < num; ii++) {
    if (inputArray[ii] > currentMax) {
      currentMax = inputArray[ii];
      currentMax_index = ii;
    }
  }

  return currentMax_index;
}