#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>

#include "output.hpp"
#include "params.hpp"

void makePlots(Params * p);

void plotPopulations(char * fileName, Params * p);

void plotKProbs(char * fileName, Params * p);

void plotCProbs(Params * p);

void plotMuFromPops(char * fileName, Params * p);

void plotDMt_z(char * fileName, Params * p);

void plotKProbsMovie(char * fileName, Params * p);

void plotCProbsMovie(char * fileName, Params * p);