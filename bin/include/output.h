#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <fftw3.h>

#include "fftmanip.h"

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputCVector(const char * fileName, fftw_complex * vec, double * x, int len);

/* Wrapper to outputCVector, which fftshifts the output.
 */
void outputCVectorShift(const char * fileName, fftw_complex * vec, double * x, int len);

#endif
