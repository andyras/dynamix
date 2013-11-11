#include "output.hpp"
#include "conversions.hpp"

//#define DEBUG_OUTPUT
#define DEBUG_OUTPUTTXPROB

/* returns true if map contains key, otherwise false */
bool isOutput(std::map<const std::string, bool> &myMap, const std::string myStr) {
  try {
    // key exists, return its value
    return myMap.at(myStr);
  }
  catch(const std::out_of_range& oor) {
    // key does not exist
    return false;
  }
}

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len) {
  // make output file
  FILE * output = NULL;
  output = fopen(fileName, "w");

  // write to the output
  for (int i = 0; i < len; i++) {
    fprintf(output, "%-.7g %-.7g %-.7g\n", x[i], vec[i][0], vec[i][1]);
  }

  // clean up the mess
  fclose(output);

  return;
}

/* Wrapper to outputFFTWVector, which fftshifts the output.
*/
void outputFFTWVectorShift(const char * fileName, fftw_complex * vec, double * x, int len) {
  // make a shifted copy of the vector to be printed
  fftw_complex * vec_shift = NULL;
  vec_shift = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*len);
  fftshift(vec, vec_shift, len);

  /*
  // make a shifted copy of the x array
  double * x_shift = new double [len];
  fftshift_double(x, x_shift, len);
  */

  // make the output
  outputFFTWVector(fileName, vec_shift, x, len);

  // clean up the mess
  fftw_free(vec_shift);
  /*
     delete [] x_shift;
     */

  return;
}

/* prints out initial wave function.  Inputs are the wave function array and
 * the number of equations.
 */
void outputWavefunction(realtype * psi, int n) {
  FILE * psi_start = NULL;
  psi_start = fopen("psi_start.out", "w");
  for (int i = 0; i < n; i++) {
    fprintf(psi_start, "%-.9g %-.9g\n", psi[i], psi[i+n]);
  }
  fclose(psi_start);
}

/* prints a vector W of length N */
void outputVector(realtype * W, int N, char * fileName) {
  FILE * out = NULL;   // output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e\n", W[ii]);
  }

  fclose(out);
}

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e %-.9e\n", evals[ii], (pow(v[ii].re,2) + pow(v[ii].im,2)));
  }

  fclose(out); 
}

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e %-.9e\n", v[ii].re, v[ii].im);
  }

  fclose(out);
}

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int jj = 0; jj < M; jj++) {
    for (int ii = 0; ii < N; ii++) {
      fprintf(out, "%-.9e %-.9e\n", v[jj*N+ii].re, v[jj*N+ii].im);
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName) {
  FILE * out = NULL; // output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    for (int jj = 0; jj < N; jj++) {
      if (jj == 0) {
	fprintf(out, "%-.9e", M[ii*N + jj]);
      }
      else {
	fprintf(out, " %-.9e", M[ii*N + jj]);
      }
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int ii = 0; ii < N; ii++) {
    fprintf(out, "%-.9e", M[ii][0]);
    for (int jj = 1; jj < N; jj++) {
      fprintf(out, " %-.9e", M[ii][jj]);
    }
    fprintf(out, "\n");
  }

  fclose(out);

  return;
}

/* Output the population in each state over time.  This function takes
 * the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputXProbsWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(fileName);

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii];
    for (int jj = start; jj < end; jj++) {
      output << " " << std::setw(8) << std::scientific
	<< (pow(wfnt[ii*p->NEQ*2 + jj],2) + pow(wfnt[ii*p->NEQ*2 + p->NEQ + jj],2));
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Output the population in each state over time.  This function takes
 * the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputXProbs(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(fileName);

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii];
    for (int jj = start; jj < end; jj++) {
      output << " "
	<< std::setw(8) << std::scientific << dmt[ii*p->NEQ2*2 + jj*p->NEQ + jj];
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Output the total population in a set of states over time */
void outputtXprobWfn(char * fileName, int start, int end, realtype * wfnt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUTTXPROB
  std::cout << "Making file " << fileName << "..." << std::endl;
  std::cout << "start index is " << start << std::endl;
  std::cout << "end index is " << end << std::endl;
#endif
  int N = p->NEQ;

  std::ofstream output(fileName);

  realtype summ;
  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ += pow(wfnt[ii*2*N + jj],2) + pow(wfnt[ii*2*N + jj + N],2);
    }
    output << std::setw(8) << std::scientific << summ << std::endl;
  }

  output.close();

#ifdef DEBUG_OUTPUTTXPROB
  std::cout << "Done making file " << fileName << "..." << std::endl;
#endif
  return;
}

/* Output the total population in a set of states over time.  This function
 * takes the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputtXprob(char * fileName, int start, int end, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "Creating file " << fileName << ".\n";
#endif
  std::ofstream output(fileName);
  realtype summ;

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << " ";
    summ = 0.0;
    for (int jj = start; jj < end; jj++) {
      summ += dmt[ii*p->NEQ2*2 + jj*p->NEQ + jj];
    }
    output << std::setw(8) << std::scientific << summ << "\n";
  }

  output.close();

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the norm of each component of the density matrix */
void outputDMZ(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting |\\rho| in time\n\n\n");
#endif

  std::ofstream output(fileName);
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj],2)
	    + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2],2));
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk],2)
	      + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2],2));
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the real part of each component of the density matrix */
void outputDMRe(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting Re(\\rho) in time\n\n\n");
#endif

  std::ofstream output(fileName);
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< dmt[2*p->NEQ2*ii + p->NEQ*jj];
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << dmt[2*p->NEQ2*ii + p->NEQ*jj + kk];
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs the imaginary part of each component of the density matrix */
void outputDMIm(char * fileName, realtype * dmt, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  fprintf(stderr, "\n\n\noutputting Im(\\rho) in time\n\n\n");
#endif

  std::ofstream output(fileName);
  // loop over time steps
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    // loop over first index
    for (int jj = 0; jj < p->NEQ; jj++) {
      // first element in row
      output << std::setw(8) << std::scientific
	<< dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2];
      // loop over second index
      for (int kk = 1; kk < p->NEQ; kk++) {
	output << " " << std::setw(8) << std::scientific
	  << dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2];
      }
      output << "\n";
    }
    output << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs energies of all states */
void outputEnergy(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(fileName);
  for (int ii = 0; ii < p->NEQ; ii++) {
    output << std::setw(8) << std::scientific << p->energies[ii] << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs all the time steps */
void outputTimes(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(fileName);
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific << p->times[ii] << "\n";
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs expectation value of energy at all times */
void outputEnergyExp(char * fileName, realtype * dmt,
    struct PARAMETERS * p) {
#ifdef DEBUG_OUTPUT
  std::cout << "\nMaking " << fileName << "\n";
#endif

  std::ofstream output(fileName);
  realtype summ ;

  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    summ = 0.0;
    for (int ii = 0; ii < p->NEQ; ii++) {
      summ += p->H[ii*p->NEQ + ii]*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii];
      for (int jj = 0; jj < ii; jj++) {
	summ += 2*p->H[ii*p->NEQ + jj]*dmt[kk*p->NEQ2*2 + jj*p->NEQ + ii];
      }
    }
    output << std::setw(8) << std::scientific
      << p->times[kk] << " "
      << std::setw(8) << std::scientific
      << summ << std::endl;
  }

#ifdef DEBUG_OUTPUT
  std::cout << "\nDone making " << fileName << std::endl;
#endif

  return;
}

/* Outputs torsional potential at simulation time points */
void outputTorsion(struct PARAMETERS * p, char * fileName) {
  std::ofstream output(fileName);

  for (int ii = 0; ii <= p->numOutputSteps; ii++) {
    output << std::setw(8) << std::scientific
      << p->times[ii] << " "
      << std::setw(8) << std::scientific
      << p->torsionV->value(p->times[ii]) << std::endl;
  }
  return;
}

/* Outputs quantities from RTA */
void outputRTA(std::map<const std::string, bool> &outs,
    realtype * dmt, struct PARAMETERS * p) {

  double ne = 0.0;
  double ekin = 0.0;
  double factor = 1.0/(M_PI*M_PI*pow(p->X2,3));
  // variables for Simpson's rule integration
  double SF = 4.0;	// Simpson's factor
  int sign = -1;	// sign
  // variables for calculating mu_e and beta
  int iter = 0;
  const int maxiter = 60;
  const double tol = 1e-12;
  const double K1 = 4.8966851;		// constants
  const double K2 = 0.04496457;
  const double K3 = 0.133376;
  const double X = 4*ne*pow(M_PI/(2*p->me),1.5)*6.9608/6.95369; // FIXME conversion at end to match Sai's values...
  double bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
  double bm = 0.0;
  const double vol = pow(1.0/5.29e-11,3);		// volume element, multiply to go from a0^-3 to m^-3
  double f = 0.0;		// value of function (f)
  double fp = 0.0;		// value of function derivative (f')
  double mue = 0.0;
  double nue = 4*ne*pow(M_PI*bm/(2*p->me),1.5);	// constant to simplify

  // vectors for time-dependent quantities
  std::vector<double> mu_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> temp_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> ne_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> ekin_t (p->numOutputSteps + 1, 0.0);
  std::vector<double> fdd_t ((p->numOutputSteps + 1)*p->Nk, 0.0);

  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    //// calculate n_e and e_kin

    // assign vector of energies
    std::vector<double> E (p->Nk,0.0);
    for (int ii = 0; ii < p->Nk; ii++) {
      E[ii] = pow(ii,2)/(2*p->me*pow(p->X2,2));
    }

    // Simpson's Rule method
    SF = 4.0;
    sign = -1;
    // skip the first point because the value will be zero
    for (int ii = 1; ii < (p->Nk-1); ii++) {
      ne += SF*factor*ii*ii*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii];
      ekin += SF*factor*pow(ii,2)*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii]*E[ii];
      SF += sign*2.0;
      sign *= -1;
    }
    // add last point
    ne += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1];
    ekin += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1]*E[p->Nk - 1];
    // divide by three
    ne /= 3.0;
    ekin /= 3.0;

    ne_t[kk] = ne;
    ekin_t[kk] = ekin;

    //// find the inverse temperature (beta)
    bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
    bm = 0.0;

    // loop applies Newton-Raphson method to get zero of function
    f = 0.0;		// value of function (f)
    fp = 0.0;		// value of function derivative (f')
    while ((fabs(bn - bm)/bm > tol) && (iter < maxiter)) {
      bm = bn;
      f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
      fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*X*pow(bm,0.5));
      bn = bm - f/fp;
      iter++;
    }

    temp_t[kk] = 1.0/bn;
    //// use beta to find chemical potential
    nue = 4*ne*pow(M_PI*bn/(2*p->me),1.5);	// constant to simplify
    mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bn;
#ifdef DEBUG_OUTPUT
    std::cout << "nue is " << nue << std::endl;
    std::cout << "(log(" << nue << ") + " << K1 << "*log(" << K2 << "*" << nue << " + 1) + " << K3 << "*" << nue << ")/" << bn << std::endl;
    std::cout << "mue is " << mue << std::endl;
#endif
    mu_t[kk] = mue;

    // update FDD
    for (int ii = 0; ii < p->Nk; ii++) {
      fdd_t[kk*p->Nk + ii] = 1.0/(1.0 + exp((E[ii] - mue)*bn));
    }
  }
  try {
    if (outs.at("mu.out")) {
      std::ofstream mu_out("mu.out");
    std::cout << "WHOOO " << std::endl << std::endl;
      for (int kk = 0; kk <= p->numOutputSteps; kk++) {
	mu_out << p->times[kk] << " " << mu_t[kk] << std::endl;
      }
      mu_out.close();

#ifdef DEBUG_OUTPUT
      std::cout << "\nDone making mu.out" << std::endl;
#endif
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    if (outs.at("temp.out")) {
      std::ofstream temp_out("temp.out");
      for (int kk = 0; kk <= p->numOutputSteps; kk++) {
	temp_out << p->times[kk] << " " << temp_t[kk]*AU2K << std::endl;
      }
      temp_out.close();

#ifdef DEBUG_OUTPUT
      std::cout << "\nDone making temp.out" << std::endl;
#endif
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    if (outs.at("ne.out")) {
      std::ofstream ne_out("ne.out");
      for (int kk = 0; kk <= p->numOutputSteps; kk++) {
	ne_out << p->times[kk] << " " << ne_t[kk] << std::endl;
      }
      ne_out.close();

#ifdef DEBUG_OUTPUT
      std::cout << "\nDone making ne.out" << std::endl;
#endif
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    if (outs.at("ekin.out")) {
      std::ofstream ekin_out("ekin.out");
      for (int kk = 0; kk <= p->numOutputSteps; kk++) {
	ekin_out << p->times[kk] << " " << ekin_t[kk] << std::endl;
      }
      ekin_out.close();

#ifdef DEBUG_OUTPUT
      std::cout << "\nDone making ekin.out" << std::endl;
#endif
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    if (outs.at("fdd.out")) {
      std::ofstream fdd_out("fdd.out");
      for (int kk = 0; kk <= p->numOutputSteps; kk++) {
	fdd_out << p->times[kk];
	for (int ii = 0; ii < p->Nk; ii++) {
	  fdd_out << " " << fdd_t[kk*p->Nk + ii];
	}
	fdd_out << std::endl;
      }
      fdd_out.close();

#ifdef DEBUG_OUTPUT
      std::cout << "\nDone making fdd.out" << std::endl;
#endif
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  return;
}

/* Outputs the laser pump intensity over time */
  void outputPumpIntensity(struct PARAMETERS * p, char * fileName) {
    std::ofstream o(fileName);
    for (int ii = 0; ii <= p->numOutputSteps; ii++) {
      o << p->times[ii] << " "
	<< gaussPulse(p->times[ii], p->pumpFWHM, p->pumpAmpl, p->pumpPeak, p->pumpFreq, p->pumpPhase)
       	<< std::endl;
    }
    o.close();

    return;
  }

/* output couplings as a matrix */
void outputCouplings(struct PARAMETERS * p, char * fileName) {
  std::ofstream o(fileName);
  //// This output is the same as the Hamiltonian, but with diagonal elements zero.
  
  // first element of first row will be zero
  o << 0;
  // rest of first row
  for (int jj = 1; jj < p->NEQ; jj++) {
    o << " " << p->H[jj];
  }
  o << std::endl;

  // remaining rows
  for (int ii = 1; ii < p->NEQ; ii++) {
    // first element in row
    o << p->H[ii*p->NEQ];
    // rest of row
    for (int jj = 1; jj < p->NEQ; jj++) {
      if (ii == jj) {
	o << " " << 0;
      }
      else {
	o << " " << p->H[ii*p->NEQ + jj];
      }
    }
    o << std::endl;
  }
  o.close();

  return;
}

/* Computes outputs independent of DM or wavefunction propagation*/
void computeGeneralOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {
  // torsion-mediated coupling
  try {
    if ((p->torsion) && (outs.at("torsion.out"))) {
      outputTorsion(p, "torsion.out");
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  
  // hamiltonian at time zero
  try {
    if (outs.at("ham.out")) {
      // &(p->H)[0] is address of first element of array in vector
      outputSquareMatrix(&(p->H)[0], p->NEQ, "ham.out");
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  
  // pump intensity
  try {
    if (outs.at("pump_intensity.out")) {
      outputPumpIntensity(p, "pump_intensity.out");
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // energies of all states
    if (outs.at("energies.out")) {
      outputEnergy("energies.out", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    if (outs.at("couplings.out")) {
      outputCouplings(p, "couplings.out");
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  return;
}

/* Computes outputs from \psi(t) */
void computeWfnOutput(realtype * wfnt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {
  // total population on all sites
  if (isOutput(outs, "totprob.out")) {
    outputtXprobWfn("totprob.out", 0, p->NEQ, wfnt, p);
  }

  // total population on bulk conduction band
  if (isOutput(outs, "tkprob.out")) {
    outputtXprobWfn("tkprob.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }

  // total population on QD
  if (isOutput(outs, "tcprob.out")) {
    outputtXprobWfn("tcprob.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }

  // total population on bridge
  if (isOutput(outs, "tbprob.out")) {
    outputtXprobWfn("tbprob.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }

  // total population on bulk valence band
  if (isOutput(outs, "tlprob.out")) {
    outputtXprobWfn("tlprob.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  // populations in all states
  if (isOutput(outs, "allprobs.out")) {
    outputXProbsWfn("allprobs.out", 0, p->NEQ, wfnt, p);
  }

  // populations in bulk CB
  if (isOutput(outs, "kprobs.out")) {
    outputXProbsWfn("kprobs.out", p->Ik, p->Ik + p->Nk, wfnt, p);
  }

  // populations in QD
  if (isOutput(outs, "cprobs.out")) {
    outputXProbsWfn("cprobs.out", p->Ic, p->Ic + p->Nc, wfnt, p);
  }

  // populations in bridge states
  if (isOutput(outs, "bprobs.out")) {
    outputXProbsWfn("bprobs.out", p->Ib, p->Ib + p->Nb, wfnt, p);
  }

  // populations in bulk VB
  if (isOutput(outs, "lprobs.out")) {
    outputXProbsWfn("lprobs.out", p->Il, p->Il + p->Nl, wfnt, p);
  }

  std::cerr << "whooooot" << std::endl;
  return;
}

/* Computes outputs from \rho(t) */
void computeDMOutput(realtype * dmt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {

  try {
    // total population
    if (outs.at("totprob.out")) {
      outputtXprob("totprob.out", 0, p->NEQ, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // populations in k states
    if (outs.at("kprobs.out")) {
      outputXProbs("kprobs.out", p->Ik, p->Ik + p->Nk, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  try {
    if (outs.at("tkprob.out")) {
      outputtXprob("tkprob.out", p->Ik, p->Ik + p->Nk, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // populations in c states
    if (outs.at("cprobs.out")) {
      outputXProbs("cprobs.out", p->Ic, p->Ic + p->Nc, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  try {
    if (outs.at("tcprob.out")) {
      outputtXprob("tcprob.out", p->Ic, p->Ic + p->Nc, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // populations in b states
    if (outs.at("bprobs.out")) {
      outputXProbs("bprobs.out", p->Ib, p->Ib + p->Nb, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  try {
    if (outs.at("tbprob.out")) {
      outputtXprob("tbprob.out", p->Ib, p->Ib + p->Nb, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // populations in l states
    if (outs.at("lprobs.out")) {
      outputXProbs("lprobs.out", p->Il, p->Il + p->Nl, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }
  try {
    if (outs.at("tlprob.out")) {
      outputtXprob("tlprob.out", p->Il, p->Il + p->Nl, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_z.out")) {
      outputDMZ("dmt_z.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_re.out")) {
      outputDMRe("dmt_re.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_im.out")) {
      outputDMIm("dmt_im.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // all time steps
    if (outs.at("times.out")) {
      outputTimes("times.out", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // expectation value of energy
    if (outs.at("energyexp.out")) {
      outputEnergyExp("energyexp.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  // RTA outputs are tied together
  outputRTA(outs, dmt, p);

  return;
}
