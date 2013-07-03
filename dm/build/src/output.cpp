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

/* Makes a gnuplot file to plot the QD populations over time */
void plot_cprobs(PARAMETERS p) {
 std::ofstream output("cprobs.plt");
 output << "#!/usr/bin/env gnuplot\n\n"
 << "reset\n"
 << "set terminal pdfcairo enhanced size 4in,3in font 'Arial-Bold,14'\n"
 << "set output '/dev/null'\n"
 << "!transpose -o _transpose ../outs/cprobs.out\n"
 << "plot '../outs/cprobs_transpose.out' every :::1 u ($1*" << p.tout << "/" << p.numOutputSteps << "):(-$2):3 matrix with image\n"
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

/* Output the population in each state over time.  This function takes
 * the indices 'start' and 'end', e.g. Ik and Ik+Nk
 */
void outputXProbs(char * fileName, int start, int end, realtype * dmt,
                  struct PARAMETERS * p) {
 std::ofstream output(fileName);

 for (int ii = 0; ii < p->numOutputSteps; ii++) {
  output << std::setw(8) << std::scientific << p->times[ii];
  for (int jj = start; jj < end; jj++) {
   output << " "
          << std::setw(8) << std::scientific << dmt[ii*p->NEQ2*2 + jj*p->NEQ + jj]
          << "\n";
  }
 }

 return;
}

/* Computes outputs from \rho(t) */
void computeDMOutput(realtype * dmt, realtype ** V, realtype * energies, realtype * t, int numTimeSteps,
                     std::map<std::string, bool> &outs, struct PARAMETERS * p) {
 // accumulator
 realtype summ;

 //// Population over time
 FILE * totprob;
 if (outs["totprob.out"]) {
  totprob = fopen("totprob.out", "w");
  for (int ii = 0; ii < numTimeSteps; ii++) {
   summ = 0.0;
   for (int jj = 0; jj < p->NEQ; jj++) {
//fprintf(stdout, "Population at time %d in state %d\n", ii, jj);
    summ += dmt[2*p->NEQ*p->NEQ*ii + p->NEQ*jj + jj];
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
   for (int jj = 0; jj < p->Nk; jj++) {
    summ += dmt[2*p->NEQ*p->NEQ*ii + p->NEQ*(p->Ik + jj) + p->Ik + jj];
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
   for (int jj = 0; jj < p->Nb; jj++) {
    summ += dmt[2*p->NEQ*p->NEQ*ii + p->NEQ*(p->Ib + jj) + p->Ib + jj];
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
   for (int jj = 0; jj < p->Nc; jj++) {
    summ += dmt[2*p->NEQ*p->NEQ*ii + p->NEQ*(p->Ic + jj) + p->Ic + jj];
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
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p->NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj],2)
                               + pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2],2)));
    for (int kk = 1; kk < p->NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", sqrt(pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk],2)
				+ pow(dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2],2)));
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
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p->NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", dmt[2*p->NEQ2*ii + p->NEQ*jj]);
    for (int kk = 1; kk < p->NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", dmt[2*p->NEQ2*ii + p->NEQ*jj + kk]);
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
  for (int ii = 0; ii < p->numOutputSteps; ii++) {
   // loop over first index
   for (int jj = 0; jj < p->NEQ; jj++) {
    // first element in row
    // loop over second index
    fprintf(dmt_z, "%+.7e", dmt[2*p->NEQ2*ii + p->NEQ*jj + p->NEQ2]);
    for (int kk = 1; kk < p->NEQ; kk++) {
     fprintf(dmt_z, " %+.7e", dmt[2*p->NEQ2*ii + p->NEQ*jj + kk + p->NEQ2]);
    }
    fprintf(dmt_z, "\n");
   }
   fprintf(dmt_z, "\n");
  }
  fclose(dmt_im);
 }

 return;
}
