#include "output.h"
#include <iostream>

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputCVector(const char * fileName, fftw_complex * vec, double * x, int len) {
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

/* Wrapper to outputCVector, which fftshifts the output.
 */
void outputCVectorShift(const char * fileName, fftw_complex * vec, double * x, int len) {
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
 outputCVector(fileName, vec_shift, x, len);

 // clean up the mess
 fftw_free(vec_shift);
 /*
 delete [] x_shift;
 */

 return;
}
