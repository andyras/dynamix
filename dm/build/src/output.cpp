#include "output.h"

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len) {
 // make output file
 FILE * output;
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
 fftw_complex * vec_shift;
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
 FILE * psi_start;
 psi_start = fopen("psi_start.out", "w");
 for (int i = 0; i < n; i++) {
  fprintf(psi_start, "%-.9g %-.9g\n", psi[i], psi[i+n]);
 }
 fclose(psi_start);
}

/* prints a vector W of length N */
void outputVector(realtype * W, int N, char * fileName) {
 int i;        // counter
 FILE * out;   // output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e\n", W[i]);
 }

 fclose(out);
}

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName) {
 int i;		// counter!
 FILE * out;	// output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e %-.9e\n", evals[i], (pow(v[i].re,2) + pow(v[i].im,2)));
 }
 
 fclose(out); 
}

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName) {
 int i;		// counter!
 FILE * out;	// output file

 out = fopen(fileName, "w");

 for (i = 0; i < N; i++) {
  fprintf(out, "%-.9e %-.9e\n", v[i].re, v[i].im);
 }

 fclose(out);
}

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName) {
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

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName) {
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

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName) {
 FILE * out;	// output file

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

/* assign coupling constants to global array V */
void buildCoupling (realtype ** vArray, PARAMETERS p,
                    std::map<std::string, bool> &outs) {
 
 int i, j;	// counters
 double Vkc;	// coupling between bulk and QD
 double Vkb1;	// coupling between bulk and first bridge
 double VbNc;	// coupling between last bridge and QD

 // initialize the coupling array
 for (i = 0; i < p.NEQ; i++) {
  for (j = 0; j < p.NEQ; j++) {
   vArray[i][j] = 0.0;
  }
 }

 // bridge
 if (p.bridge_on) {
  // coupling between k and b1
  if ((p.scale_bubr) && (p.Nk > 1)) {
   Vkb1 = sqrt(p.Vbridge[0]*(p.kBandTop-p.kBandEdge)/(p.Nk-1));
  }
  else {
   Vkb1 = p.Vbridge[0];
  }
  if (p.parabolicCoupling) {
   for (i = 0; i < p.Nk; i++) {
    vArray[p.Ik+i][p.Ib] = parabolicV(Vkb1, p.energies[p.Ik+i], p.kBandEdge, p.kBandTop);
    vArray[p.Ib][p.Ik+i] = parabolicV(Vkb1, p.energies[p.Ik+i], p.kBandEdge, p.kBandTop);
   }
  }
  else {
   for (i = 0; i < p.Nk; i++) {
    vArray[p.Ik+i][p.Ib] = Vkb1;
    vArray[p.Ib][p.Ik+i] = Vkb1;
   }
  }
   
  // coupling between bN and c
  if ((p.scale_brqd) && (p.Nc > 1)) {
   VbNc = p.Vbridge[p.Nb]/sqrt(p.Nc-1);
  }
  else {
   VbNc = p.Vbridge[p.Nb];
  }
  for (i = 0; i < p.Nc; i++) {
   vArray[p.Ic+i][p.Ib+p.Nb-1] = VbNc;
   vArray[p.Ib+p.Nb-1][p.Ic+i] = VbNc;
  }
  
  // coupling between bridge states
  for (i = 0; i < p.Nb - 1; i++) {
   vArray[p.Ib+i][p.Ib+i+1] = p.Vbridge[i+1];
   vArray[p.Ib+i+1][p.Ib+i] = p.Vbridge[i+1];
  }
 }
 // no bridge
 else {				
  // scaling
  if ((p.scale_buqd) && (p.Nk > 1)) {
   Vkc = sqrt(p.Vnobridge[0]*(p.kBandTop-p.kBandEdge)/(p.Nk-1));
  }
  else {
   Vkc = p.Vnobridge[0];
  }

  // parabolic coupling of bulk band to QD
  if (p.parabolicCoupling) {
   for (i = 0; i < p.Nk; i++) {
    for (j = 0; j < p.Nc; j++) {
     vArray[p.Ik+i][p.Ic+j] = parabolicV(Vkc, p.energies[p.Ik+i], p.kBandEdge, p.kBandTop);
     vArray[p.Ic+j][p.Ik+i] = parabolicV(Vkc, p.energies[p.Ik+i], p.kBandEdge, p.kBandTop);
    }
   }
  }
  else {
   for (i = 0; i < p.Nk; i++) {
    for (j = 0; j < p.Nc; j++) {
     vArray[p.Ik+i][p.Ic+j] = Vkc;
     vArray[p.Ic+j][p.Ik+i] = Vkc;
    }
   }
  }
 }

#ifdef DEBUG
 cout << "\nCoupling matrix:\n";
 for (i = 0; i < p.NEQ; i++) {
  for (j = 0; j < p.NEQ; j++)
   cout << scientific << vArray[i][j] << " ";
  cout << endl;
 }
#endif

 FILE * couplings;
 if (outs["couplings.out"]) {
  couplings = fopen("couplings.out","w");
  for (i = 0; i < p.NEQ; i++) {
   for (j = 0; j < p.NEQ; j++) {
    fprintf(couplings,"%.7g ",vArray[i][j]);
   }
   fprintf(couplings,"\n");
  }
  fclose(couplings);
 }
 
}


/* Computes outputs from \rho(t) */
void computeDMOutput(realtype * dmt, realtype ** V, realtype * energies, realtype * t, int numTimeSteps,
                     std::map<std::string, bool> &outs, PARAMETERS p) {
 // accumulator
 realtype summ;

 //// Population over time
 FILE * totprob;
 if (outs["totprob.out"]) {
  totprob = fopen("totprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < p.NEQ; jj++) {
//fprintf(stdout, "Population at time %d in state %d\n", ii, jj);
    summ += dmt[2*p.NEQ*p.NEQ*ii + p.NEQ*jj + jj];
   }
   fprintf(totprob, "%-.7g %-.7g\n", t[ii], summ);
  }
  fclose(totprob);
 }

 //// Population in k states
 FILE * tkprob;
 if (outs["tkprob.out"]) {
  tkprob = fopen("tkprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < p.Nk; jj++) {
    summ += dmt[2*p.NEQ*p.NEQ*ii + p.NEQ*(p.Ik + jj) + p.Ik + jj];
   }
   fprintf(tkprob, "%-.7g %-.7g\n", t[ii], summ);
  }
  fclose(tkprob);
 }

 //// Population in b states
 FILE * tbprob;
 if (outs["tbprob.out"]) {
  tbprob = fopen("tbprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < p.Nb; jj++) {
    summ += dmt[2*p.NEQ*p.NEQ*ii + p.NEQ*(p.Ib + jj) + p.Ib + jj];
   }
   fprintf(tbprob, "%-.7g %-.7g\n", t[ii], summ);
  }
  fclose(tbprob);
 }

 //// Population in c states
 FILE * tcprob;
 if (outs["tcprob.out"]) {
  tcprob = fopen("tcprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < p.Nc; jj++) {
    summ += dmt[2*p.NEQ*p.NEQ*ii + p.NEQ*(p.Ic + jj) + p.Ic + jj];
   }
   fprintf(tcprob, "%-.7g %-.7g\n", t[ii], summ);
  }
  fclose(tcprob);
 }

#ifdef DEBUG
 fprintf(stderr, "\n\n\noutputting dm in time\n\n\n");
#endif

 // norm of DM elements
 FILE * dmt_z;
 if (outs["dmt_z.out"]) {
  dmt_z = fopen("dmt_z.out", "w");
  // loop over time steps
  for (int ii = 0; ii < p.numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p.NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", sqrt(pow(dmt[2*p.NEQ2*ii + p.NEQ*jj],2)
                               + pow(dmt[2*p.NEQ2*ii + p.NEQ*jj + p.NEQ2],2)));
    for (int kk = 1; kk < p.NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", sqrt(pow(dmt[2*p.NEQ2*ii + p.NEQ*jj + kk],2)
				+ pow(dmt[2*p.NEQ2*ii + p.NEQ*jj + kk + p.NEQ2],2)));
    }
    fprintf(dmt_z, "\n");
   }
   fprintf(dmt_z, "\n");
  }
  fclose(dmt_z);
 }

 // real part of DM
 FILE * dmt_re;
 if (outs["dmt_re.out"]) {
  dmt_re = fopen("dmt_re.out", "w");
  // loop over time steps
  for (int ii = 0; ii < p.numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p.NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", dmt[2*p.NEQ2*ii + p.NEQ*jj]);
    for (int kk = 1; kk < p.NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", dmt[2*p.NEQ2*ii + p.NEQ*jj + kk]);
    }
    fprintf(dmt_z, "\n");
   }
   fprintf(dmt_z, "\n");
  }
  fclose(dmt_re);
 }

 // imaginary part of DM
 FILE * dmt_im;
 if (outs["dmt_im.out"]) {
  dmt_im = fopen("dmt_im.out", "w");
  // loop over time steps
  for (int ii = 0; ii < p.numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p.NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", dmt[2*p.NEQ2*ii + p.NEQ*jj + p.NEQ2]);
    for (int kk = 1; kk < p.NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", dmt[2*p.NEQ2*ii + p.NEQ*jj + kk + p.NEQ2]);
    }
    fprintf(dmt_z, "\n");
   }
   fprintf(dmt_z, "\n");
  }
  fclose(dmt_im);
 }

 return;
}
