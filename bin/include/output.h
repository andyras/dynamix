#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <fftw3.h>

#include "fftmanip.h"

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len);

/* Wrapper to outputFFTWVector, which fftshifts the output.
 */
void outputFFTWVectorShift(const char * fileName, fftw_complex * vec, double * x, int len);

#endif
