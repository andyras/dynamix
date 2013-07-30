#include "output.hpp"

//#define DEBUG_OUTPUT
//#define DEBUG_OUTPUTTXPROB

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
  int i;        // counter
  FILE * out = NULL;   // output file

  out = fopen(fileName, "w");

  for (i = 0; i < N; i++) {
    fprintf(out, "%-.9e\n", W[i]);
  }

  fclose(out);
}

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName) {
  int i;		// counter!
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (i = 0; i < N; i++) {
    fprintf(out, "%-.9e %-.9e\n", evals[i], (pow(v[i].re,2) + pow(v[i].im,2)));
  }

  fclose(out); 
}

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName) {
  int i;		// counter!
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (i = 0; i < N; i++) {
    fprintf(out, "%-.9e %-.9e\n", v[i].re, v[i].im);
  }

  fclose(out);
}

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName) {
  int i, j;	// counters!
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (j = 0; j < M; j++) {
    for (i = 0; i < N; i++) {
      fprintf(out, "%-.9e %-.9e\n", v[j*N+i].re, v[j*N+i].im);
    }
    fprintf(out, "\n");
  }

  fclose(out);
}

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName) {
  int i, j;     // counters
  FILE * out = NULL; // output file

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

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName) {
  FILE * out = NULL;	// output file

  out = fopen(fileName, "w");

  for (int i = 0; i < N; i++) {
    fprintf(out, "%-.9e", M[i][0]);
    for (int j = 1; j < N; j++) {
      fprintf(out, " %-.9e", M[i][j]);
    }
    fprintf(out, "\n");
  }

  fclose(out);

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
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // populations in k states
    if (outs.at("kprobs.out")) {
      outputXProbs("kprobs.out", p->Ik, p->Ik + p->Nk, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }
  try {
    if (outs.at("tkprob.out")) {
      outputtXprob("tkprob.out", p->Ik, p->Ik + p->Nk, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // populations in c states
    if (outs.at("cprobs.out")) {
      outputXProbs("cprobs.out", p->Ic, p->Ic + p->Nc, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }
  try {
    if (outs.at("tcprob.out")) {
      outputtXprob("tcprob.out", p->Ic, p->Ic + p->Nc, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // populations in b states
    if (outs.at("bprobs.out")) {
      outputXProbs("bprobs.out", p->Ib, p->Ib + p->Nb, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }
  try {
    if (outs.at("tbprob.out")) {
      outputtXprob("tbprob.out", p->Ib, p->Ib + p->Nb, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // populations in l states
    if (outs.at("lprobs.out")) {
      outputXProbs("lprobs.out", p->Il, p->Il + p->Nl, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }
  try {
    if (outs.at("tlprob.out")) {
      outputtXprob("tlprob.out", p->Il, p->Il + p->Nl, dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_z.out")) {
      outputDMZ("dmt_z.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_re.out")) {
      outputDMRe("dmt_re.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // norm of DM elements
    if (outs.at("dmt_im.out")) {
      outputDMIm("dmt_im.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // energies of all states
    if (outs.at("energies.out")) {
      outputEnergy("energies.out", p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // all time steps
    if (outs.at("times.out")) {
      outputTimes("times.out", p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  try {
    // expectation value of energy
    if (outs.at("energyexp.out")) {
      outputEnergyExp("energyexp.out", dmt, p);
    }
  }
  catch (const std::out_of_range& oor) {
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
  }

  return;
}
