#include "output.h"
#include <iostream>

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

/* prints out the density matrix varying in time to files */
void outputDMt(realtype * dmt, int NEQ, int numOutputSteps, std::map<std::string, bool> &outs) {
fprintf(stderr, "\n\n\noutputting dm in time\n\n\n");
 FILE * dmt_z;
 if (outs["dmt_z.out"]) {
  dmt_z = fopen("dmt_z.out", "w");
  // loop over time steps
  for (int ii = 0; ii < numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", dmt[2*NEQ*NEQ*ii + NEQ*jj]);
    for (int kk = 1; kk < NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", dmt[2*NEQ*NEQ*ii + NEQ*jj+ kk]);
    }
    fprintf(dmt_z, "\n");
   }
   fprintf(dmt_z, "\n");
  }
 }
 fclose(dmt_z);

 if (outs["dmt_re.out"]) {
 }

 if (outs["dmt_im.out"]) {
 }
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

/* Computes outputs from \rho(t) */
void computeDMOutput(realtype * dmt, int NEQ, realtype ** V, realtype * energies, realtype * t, int numTimeSteps,
                     std::map<std::string, bool> &outs) {
 // accumulator
 realtype summ;

 //// Population over time
 FILE * totprob;
 if (outs["totprob.out"]) {
  totprob = fopen("totprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < NEQ; jj++) {
//fprintf(stdout, "Population at time %d in state %d\n", ii, jj);
    summ += dmt[2*NEQ*NEQ*ii + NEQ*jj + jj];
   }
   fprintf(totprob, "%-.7g %-.7g\n", t[ii], summ);
  }
  fclose(totprob);
 }

 //// Population in k states

 //// Population in b states

 //// Population in c states

 return;
}
