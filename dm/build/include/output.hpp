#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <map>
#include <cmath>
#include <string>
#include <nvector/nvector_serial.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "conversions.hpp"
#include "numerical.hpp"
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
    const realtype * dmt, Params * p);

void outputIntegratedWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, Params * p);

void outputIntegralDM(const std::string fileName, const int start, const int end,
    const realtype * dmt, Params * p);

void outputIntegralWfn(const std::string fileName, const int start, const int end,
    const realtype * wfnt, Params * p);

void outputXProbsWfn(char * fileName, int start, int end, realtype * wfnt,
    Params * p);

void outputXProbs(char * fileName, int start, int end, realtype * dmt,
    Params * p);

void outputtXprobWfn(char * fileName, int start, int end, realtype * wfnt,
    Params * p);

void outputtXprob(char * fileName, int start, int end, realtype * dmt,
    Params * p);

void outputDMZ(realtype * dmt, Params * p);

void outputDMCoherences(realtype * dmt, Params * p);

void outputDMIm(char * fileName, realtype * dmt, Params * p);

void outputDMRe(char * fileName, realtype * dmt, Params * p);

void outputEnergy(char * fileName, Params * p);

void outputEnergyExpWfn(const char * fileName, Params * p);

void outputTimes(char * fileName, Params * p);

void outputEnergyExp(char * fileName, realtype * dmt,
    Params * p);

void outputDynamicMu(char * fileName, realtype * dmt, int bandFlag, Params * p);

void outputMuFromPops(char * fileName, realtype * dmt,
    Params * p);

void outputTorsion(Params * p, char * fileName);

void outputRTA(std::map<const std::string, bool> &outs,
    realtype * dmt, Params * p);

void outputCouplings(Params * p, char * fileName);

void outputTorsionSin2(Params * p, char * fileName);

void findPeaksWfn(char * fileName, int start, int end, realtype * wfnt,
    Params * p);

void outputDeriv(char * fileName, int n, realtype * deriv, Params * p);

void outputDerivsWfn(std::map<const std::string, bool> &outs, realtype * wfnt,
    Params * p);

void outputDerivsDM(std::map<const std::string, bool> &outs, realtype * dmt,
    Params * p);

void computeGeneralOutputs(Params * p);

void computeWfnOutput(realtype * wfnt, Params * p);

void computeDMOutput(realtype * dmt, Params * p);
#endif
