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

/* Outputs torsional potential at simulation time points */
void outputTorsion(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p, char * fileName) {
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
void outputRTA(realtype * dmt, std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {

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
    //#ifdef DEBUG_RTA
    //std::cout << "XX " << X/pow(2.293710449e+17,1.5) << std::endl;
    //#endif
    double bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
    double bm = 0.0;
    const double vol = pow(1.0/5.29e-11,3);		// volume element, multiply to go from a0^-3 to m^-3
    double f = 0.0;		// value of function (f)
    double fp = 0.0;		// value of function derivative (f')
    double mue = 0.0;
    double nue = 4*ne*pow(M_PI*bm/(2*p->me),1.5);	// constant to simplify
  
  // loop over timesteps
  for (int kk = 0; kk <= p->numOutputSteps; kk++) {
    //#ifdef DEBUG_RTA
    //std::cout << "p->X2   " << p->X2 << std::endl;
    //#endif

    //// calculate n_e and e_kin
    //
    //#ifdef DEBUG_RTA
    //std::cout << "factor   " << factor << std::endl;
    //#endif
    // assign vector of energies
    std::vector<double> E (p->Nk,0.0);
    for (int ii = 0; ii < p->Nk; ii++) {
      E[ii] = pow(ii,2)/(2*p->me*pow(p->X2,2));
    }

    //#ifdef DEBUG_RTA
    //std::cout << std::setprecision(28);
    //#endif
    // Simpson's Rule method
    SF = 4.0;
    sign = -1;
    // skip the first point because the value will be zero
    for (int ii = 1; ii < (p->Nk-1); ii++) {
      //ne += SF*factor*ii*ii*NV_Ith_S(y, ii*p->NEQ + ii);
      ne += SF*factor*ii*ii*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii];
      //ekin += SF*factor*pow(ii,2)*NV_Ith_S(y, ii*p->NEQ + ii)*E[ii];
      ekin += SF*factor*pow(ii,2)*dmt[kk*p->NEQ2*2 + ii*p->NEQ + ii]*E[ii];
      //#ifdef DEBUG_RTA
      //std::cout << "Ne " << ii*ii << "*" << SF << "/3.0*" << NV_Ith_S(y, ii*p->NEQ + ii) << "/" << pow(5.29e-11,3)/factor << std::endl;
      //std::cout << "ekin " << pow(ii,4) << "*" << SF << "/3.0*" << NV_Ith_S(y, ii*p->NEQ + ii) << "*" << 4.3597482e-18/(2*p->me*pow(p->X2,2)) << "/" << pow(5.29e-11,3)/factor << std::endl;
      //std::cout << "ekin " << ekin/pow(5.29e-11,3)*4.3597482e-18/3.0 
      //<< " += " << SF*factor*pow(ii,4)*NV_Ith_S(y, ii*p->NEQ + ii)/(2*p->me*pow(p->X2,2))/pow(5.29e-11,3)*4.3597482e-18/3.0 << std::endl;
      //#endif
      SF += sign*2.0;
      sign *= -1;
    }
    // add last point
    ne += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1];
    ekin += factor*pow(p->Nk-1,2)*dmt[kk*p->NEQ2*2 + (p->Nk - 1)*p->NEQ + p->Nk - 1]*E[p->Nk - 1];
    // divide by three
    ne /= 3.0;
    ekin /= 3.0;

    //#ifdef DEBUG_RTA
    std::cout << "ne   (" << kk << ") " << ne << std::endl;
    //std::cout << "ne (SI)   " << ne/pow(5.29e-11,3) << std::endl;
    std::cout << "ekin (" << kk << ") " << ekin << std::endl;
    //std::cout << "ekin (SI) " << ekin/pow(5.29e-11,3)*4.3597482e-18 << std::endl;
    //#endif

    //// find the inverse temperature (beta)
    bn = 1.9e20*4.3597482e-18*0.5;		// intermediate values of beta; bn is higher iteration
    bm = 0.0;

    // loop applies Newton-Raphson method to get zero of function
    f = 0.0;		// value of function (f)
    fp = 0.0;		// value of function derivative (f')
    //#ifdef DEBUG_RTA
    //std::cout << "Newton-Raphson to find inverse temperature" << std::endl;
    //#endif
    while ((fabs(bn - bm)/bm > tol) && (iter < maxiter)) {
      bm = bn;
      f = -bm*ekin + 1.5*ne*(1 + K1 - K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5)) + 0.5*K3*X*pow(bm,1.5));
      fp = -ekin + 2.25*ne*(K1/(K2*X*pow(bm,2.5))*log(1 + K2*X*pow(bm,1.5)) - K1/(bm*(1 + K2*X*pow(bm,1.5))) + 0.5*K3*X*pow(bm,0.5));
      //#ifdef DEBUG_RTA
      //std::cout << "Iteration     " << std::setw(15) << iter << std::endl;
      //std::cout << "bm            " << std::setw(15) << bm/4.3597482e-18 << std::endl;
      //std::cout << "f(bm) term 1: " << std::setw(15) << vol*-bm*ekin << std::endl;
      //std::cout << "f(bm) term 2: " << std::setw(15) << vol*1.5*ne*(1 + K1) << std::endl;
      //std::cout << "f(bm) term 3: " << std::setw(15) << vol*1.5*ne*(-1*K1/(K2*X)*pow(bm,-1.5)*log(1 + K2*X*pow(bm,1.5))) << std::endl;
      //std::cout << "f(bm) term 4: " << std::setw(15) << vol*1.5*ne*(0.5*K3*X*pow(bm,1.5)) << std::endl;
      //std::cout << "f(bm) (SI)    " << std::setw(15) << f*pow(1.0/5.29e-11,3) << std::endl;
      //std::cout << "f(bm) (a.u)   " << std::setw(15) << f << std::endl;
      //std::cout << "f'(bm) (SI)   " << std::setw(15) << fp*pow(1.0/5.29e-11,3)*4.3597482e-18 << std::endl;
      //std::cout << "f'(bm) (a.u)  " << std::setw(15) << fp << std::endl;
      //#endif
      bn = bm - f/fp;
      iter++;
    }
    //#ifdef DEBUG_RTA
    //std::cout << std::endl;
    //#endif

    //// use beta to find chemical potential
    nue = 4*ne*pow(M_PI*bm/(2*p->me),1.5);	// constant to simplify
    mue = (log(nue) + K1*log(K2*nue + 1) + K3*nue)/bm;
    //#ifdef DEBUG_RTA
    //std::cout << "Chemical potential " << mue*4.3597482e-18 << std::endl;
    //#endif

    //#ifdef DEBUG_RTA
    //std::cout << "inverse temp is " << bn << std::endl;
    //std::cout << std::endl;
    //#endif
    /*
    std::ofstream fddout("fdd.out");
    for (int ii = 0; ii < p->Nk; ii++) {
      fdd[ii] = 1.0/(1.0 + exp((E[ii] - mue)*bn));
      fddout << E[ii]*27.211 << " " << fdd[ii] << std::endl;
      //#ifdef DEBUG_RTA
      //std::cout << "FDD[" << ii << "]: " << std::scientific << fdd[ii] << std::endl;
      //#endif
    }
    fddout.close();
    */
  }

  return;
}

/* Computes outputs independent of DM or wavefunction */
void computeGeneralOutputs(std::map<const std::string, bool> &outs,
    struct PARAMETERS * p) {
  try {
    if ((p->torsion) && (outs.at("torsion.out"))) {
      outputTorsion(outs, p, "torsion.out");
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG_OUTPUT
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
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
  outputRTA(dmt, outs, p);

  return;
}
