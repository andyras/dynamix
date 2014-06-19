#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <map>
#include <cmath>
#include <string>
#include <fftw3.h>
#include <nvector/nvector_serial.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "conversions.hpp"
#include "numerical.hpp"
#include "fftmanip.hpp"
#include "dynamix.hpp"
#include "params.hpp"
#include "rhs.hpp"

/* TYPE DEFINITIONS */

/* structure to hold complex numbers for linear algebra routines */
typedef struct {
  double re;
  double im;
} complex16;

/* FUNCTIONS */

bool isOutput(std::map<const std::string, bool> &myMap, const std::string myStr);

std::string outputFileName(char * fileName, Params * p);

/* prints out array of fftw_complex values.  The 'x' array is
 * the x-axis variable: time, energy, &c.
 */
void outputFFTWVector(const char * fileName, fftw_complex * vec, double * x, int len);

/* Wrapper to outputFFTWVector, which fftshifts the output.
*/
void outputFFTWVectorShift(const char * fileName, fftw_complex * vec, double * x, int len);

/* prints out initial wave function.  Inputs are the wave function array and
 * the number of equations.
 */
void outputWavefunction(realtype * psi, int n);

/* prints a vector W of length N */
void outputVector(realtype * W, int N, char * fileName);

/* prints the elementwise square of a complex vector */
void outputPsiSquare(complex16 * v, realtype * evals,  int N, char * fileName);

/* prints a complex vector v of length N */
void outputCVector(complex16 * v, int N, char * fileName);

/* prints a complex vector v of length N with M time steps */
void outputCVectorTime(complex16 * v, int N, int M, char * fileName);

/* prints a square matrix M of dimension N */
void outputSquareMatrix(realtype * M, int N, char * fileName);

/* prints a square matrix stored as a 2D array */
void output2DSquareMatrix(realtype ** M, int N, char * fileName);

void outputIntegratedDM(const std::string fileName, const int start, const int end,
    const realtype * dmt, struct Params * p);

void outputIntegratedWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, struct Params * p);

void outputIntegralDM(const std::string fileName, const int start, const int end,
    const realtype * dmt, struct Params * p);

void outputIntegralWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, struct Params * p);

void outputXProbsWfn(char * fileName, int start, int end, realtype * wfnt,
    struct Params * p);

void outputXProbs(char * fileName, int start, int end, realtype * dmt,
    struct Params * p);

void outputtXprobWfn(char * fileName, int start, int end, realtype * wfnt,
    struct Params * p);

void outputtXprob(char * fileName, int start, int end, realtype * dmt,
    struct Params * p);

void outputDMZ(realtype * dmt, struct Params * p);

void outputDMCoherences(realtype * dmt, struct Params * p);

void outputDMIm(char * fileName, realtype * dmt, struct Params * p);

void outputDMRe(char * fileName, realtype * dmt, struct Params * p);

void outputEnergy(char * fileName, struct Params * p);

void outputEnergyExpWfn(const char * fileName, struct Params * p);

void outputTimes(char * fileName, struct Params * p);

void outputEnergyExp(char * fileName, realtype * dmt,
    struct Params * p);

void outputDynamicMu(char * fileName, realtype * dmt, int bandFlag, struct Params * p);

void outputMuFromPops(char * fileName, realtype * dmt,
    struct Params * p);

void outputTorsion(struct Params * p, char * fileName);

void outputRTA(std::map<const std::string, bool> &outs,
    realtype * dmt, struct Params * p);

void outputCouplings(struct Params * p, char * fileName);

void outputTorsionSin2(struct Params * p, char * fileName);

void findPeaksWfn(char * fileName, int start, int end, realtype * wfnt,
    struct Params * p);

void outputDeriv(char * fileName, int n, realtype * deriv, struct Params * p);

void outputDerivsWfn(std::map<const std::string, bool> &outs, realtype * wfnt,
    struct Params * p);

void outputDerivsDM(std::map<const std::string, bool> &outs, realtype * dmt,
    struct Params * p);

void computeGeneralOutputs(std::map<const std::string, bool> &outs,
    struct Params * p);

void computeWfnOutput(realtype * wfnt, std::map<const std::string, bool> &outs,
    struct Params * p);

void computeDMOutput(realtype * dmt, std::map<const std::string, bool> &outs,
    struct Params * p);
#endif
