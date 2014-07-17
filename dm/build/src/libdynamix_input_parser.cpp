#include "libdynamix_input_parser.hpp"

// #define DEBUG
// #define DEBUG_INPUT_PARSER

void assignOutputs(const char * inputFile, std::map<const std::string, bool> &outs,
    Params * p) {
  std::string line;
  std::ifstream input;

  input.open(inputFile, std::ios::in);

#ifdef DEBUG_INPUT_PARSER
  std::cout << "Parsing file " << inputFile << "\n";
#endif

  if (input.good() == false) {
    std::cerr << "[assignOutputs] ERROR: bad input file " << inputFile << "\n";
  }

  // get first line
  getline(input, line);

  // skip lines until header specifying start of outputs
  while (getline(input, line)) {
#ifdef DEBUG_INPUT_PARSER
    std::cout << "line is " << line << "\n";
#endif
    if (line != "[[Output Files]]") {
      continue;
    }
    else {
      break;
    }
  }

  // read inputs until [[End]] line or EOF
  while (getline(input, line)) {
    if (line != "[[End]]") {

      // skip comments
      if (line.substr(0,1) == "#") {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Skipping comment line: " << line << "\n";
#endif
	continue;
      }

      // skip whitespace (space/tab) and blank lines
      else if ((line.find_first_not_of(' ') == std::string::npos)
	  || (line.find_first_not_of('\t') == std::string::npos)) {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Skipping blank/whitespace line\n";
#endif
	continue;
      }

      // turn on outputs which are in this section of the input
      else {
#ifdef DEBUG_INPUT_PARSER
	std::cout << "Creating output file:  " << line << "\n";
#endif
	outs.insert(std::pair<const std::string,bool>(line,true));
	if ((line.substr(line.length()-4, line.length()) != ".out")
	    && (line.substr(line.length()-4, line.length()) != ".plt")) {
	  std::cerr << "WARNING [" << __FUNCTION__ << "]: output file extension is not '.out' or '.plt'\n";
	}
      }
    }
  }

  // if creating plot files, make sure that the appropriate outputs are being created.
  if (isOutput(outs, "projections.plt")) {
    outs.insert(std::pair<const std::string,bool>("tkprob.out", true));
    outs.insert(std::pair<const std::string,bool>("tcprob.out", true));

    // create bridge output if bridge is on
    if (p->bridge_on) {
      outs.insert(std::pair<const std::string,bool>("tbprob.out", true));
    }
  }

  if (isOutput(outs, "cprobs.plt")) {
    outs.insert(std::pair<const std::string,bool>("cprobs.out", true));
  }

  if ((isOutput(outs, "kprobs.plt")) || (isOutput(outs, "kprobs_movie.plt"))) {
    outs.insert(std::pair<const std::string,bool>("kprobs.out", true));
  }

  if ((isOutput(outs, "dmt_z.plt")) && (! p->wavefunction)) {
    outs.insert(std::pair<const std::string,bool>("dmt_z.out", true));
  }

  // create output directory if it does not exist //////////////////////////////

  if (!fileExists(p->outputDir)) {
    mkdir(p->outputDir.c_str(), 0755);
  }
}

/* assigns params to the Params struct from the input file */
void assignParams(std::string inputFile, Params * p) {
  std::string line;
  std::string input_param;
  std::string param_val;
  size_t equals_pos;
  size_t space_pos;

#ifdef DEBUG_INPUT_PARSER
  std::cout << "\nParsing file: " << inputFile << "\n";
#endif

  std::ifstream bash_in;  // declare input file stream

  bash_in.open(inputFile.c_str(), std::ios::in);	// open file as input stream
  if (bash_in.good() == false) {
    fprintf(stderr, "ERROR [Inputs]: file '%s' not available for reading\n", inputFile.c_str());
    exit(-1);
  }

  // read first line of input file
  getline (bash_in,line);

  // skip non-parameter lines
  while ( line != "## START INPUT PARAMETERS ##") {
    getline (bash_in,line);
  }

  while ( line != "## END INPUT PARAMETERS ##") {
    // skip comment lines
    if ( line.substr(0,1) == "#" ) {
      getline (bash_in,line);
      continue;
    }
    // find first equals sign
    equals_pos=line.find("=");
    // find first whitespace
    space_pos=(line.find(" ") > line.find("\t") ? line.find("\t") : line.find(" "));
    // parameter name is before equals sign
    input_param = line.substr(0,int(equals_pos));
    // parameter is after equals sign, before space
    param_val = line.substr(int(equals_pos)+1,int(space_pos)-int(equals_pos));
    // extract parameters
#ifdef DEBUG
    std::cout << "Parameter: " << input_param << std::endl << "New value: " << stof(param_val.c_str()) << std::endl;
#endif
    if (input_param == "timedepH") { p->timedepH = stoi(param_val); }
    else if (input_param == "justPlots") { p->justPlots = stoi(param_val); }
    else if (input_param == "nproc") { p->nproc = stoi(param_val); }
    else if (input_param == "wavefunction") { p->wavefunction = stoi(param_val); }
    else if (input_param == "coherent") { p->coherent = stoi(param_val); }
    else if (input_param == "analytical") { p->analytical = stoi(param_val); }
    else if (input_param == "kinetic") { p->kinetic = stoi(param_val); }
    else if (input_param == "kineticQD") { p->kineticQD = stoi(param_val); }
    else if (input_param == "dynamicMu") { p->dynamicMu = stoi(param_val); }
    else if (input_param == "dynamicMuQD") { p->dynamicMuQD = stoi(param_val); }
    else if (input_param == "rta") { p->rta = stoi(param_val); }
    else if (input_param == "rtaQD") { p->rtaQD = stoi(param_val); }
    else if (input_param == "dephasing") { p->dephasing = stoi(param_val); }
    else if (input_param == "progressFile") { p->progressFile = stoi(param_val); }
    else if (input_param == "progressStdout") { p->progressStdout = stoi(param_val); }
    else if (input_param == "abstol") { p->abstol = stof(param_val); }
    else if (input_param == "reltol" ) { p->reltol = stof(param_val); }
    else if (input_param == "tout" ) { p->tout = stof(param_val); }
    else if (input_param == "numsteps" ) { p->numsteps = stoi(param_val); }
    else if (input_param == "numOutputSteps" ) { p->numOutputSteps = stoi(param_val); }
    else if (input_param == "kBandEdge" ) { p->kBandEdge = stof(param_val); }
    else if (input_param == "kBandTop" ) { p->kBandTop = stof(param_val); }
    else if (input_param == "lBandTop" ) { p->lBandTop = stof(param_val); }
    else if (input_param == "bulk_gap" ) { p->bulk_gap = stof(param_val); }
    else if (input_param == "Nk" ) { p->Nk = stoi(param_val); }
    else if (input_param == "Nk_first" ) { p->Nk_first = stoi(param_val); }
    else if (input_param == "Nk_final" ) { p->Nk_final = stoi(param_val); }
    else if (input_param == "Nc_first" ) { p->Nc_first = stoi(param_val); }
    else if (input_param == "Nc_final" ) { p->Nc_final = stoi(param_val); }
    else if (input_param == "valenceBand" ) { p->valenceBand = stof(param_val); }
    else if (input_param == "Nl" ) { p->Nl = stoi(param_val); }
    else if (input_param == "bulkGaussSigma" ) { p->bulkGaussSigma = stof(param_val); }
    else if (input_param == "bulkGaussMu" ) { p->bulkGaussMu = stof(param_val); }
    else if (input_param == "me" ) { p->me = stof(param_val); }
    else if (input_param == "mh" ) { p->mh = stof(param_val); }
    else if (input_param == "me_c" ) { p->me_c = stof(param_val); }
    else if (input_param == "mh_c" ) { p->mh_c = stof(param_val); }
    else if (input_param == "X2" ) { p->X2 = stof(param_val); }
    else if (input_param == "temperature" ) { p->temperature = stof(param_val); }
    else if (input_param == "EF" ) { p->EF = stof(param_val); }
    else if (input_param == "gamma1" ) { p->gamma1 = stof(param_val); }
    else if (input_param == "gamma1_c" ) { p->gamma1_c = stof(param_val); }
    else if (input_param == "gamma2" ) { p->gamma2 = stof(param_val); }
    else if (input_param == "muLK" ) { p->muLK = stof(param_val); }
    else if (input_param == "pumpFWHM" ) { p->pumpFWHM = stof(param_val); }
    else if (input_param == "pumpPeak" ) { p->pumpPeak = stof(param_val); }
    else if (input_param == "pumpFreq" ) { p->pumpFreq = stof(param_val); }
    else if (input_param == "pumpAmpl" ) { p->pumpAmpl = stof(param_val); }
    else if (input_param == "pumpPhase" ) { p->pumpPhase = stof(param_val); }
    else if (input_param == "CBPopFlag" ) { p->CBPopFlag = stoi(param_val); }
    else if (input_param == "VBPopFlag" ) { p->VBPopFlag = stoi(param_val); }
    else if (input_param == "QDPopFlag" ) { p->QDPopFlag = stoi(param_val); }
    else if (input_param == "laser_on" ) { p->laser_on = stoi(param_val); }
    else if (input_param == "parabolicCoupling" ) { p->parabolicCoupling = stoi(param_val); }
    else if (input_param == "scale_bubr" ) { p->scale_bubr = stoi(param_val); }
    else if (input_param == "scale_brqd" ) { p->scale_brqd = stoi(param_val); }
    else if (input_param == "scale_buqd" ) { p->scale_buqd = stoi(param_val); }
    else if (input_param == "scale_laser" ) { p->scale_laser = stoi(param_val); }
    else if (input_param == "bridge_on" ) { p->bridge_on = stoi(param_val); }
    else if (input_param == "random_phase" ) { p->random_phase = stoi(param_val); }
    else if (input_param == "random_seed" ) { p->random_seed = stoi(param_val); }
    else if (input_param == "torsion" ) { p->torsion = stoi(param_val); }
    else if (input_param == "torsionFile" ) { p->torsionFile = param_val; }
    else if (input_param == "torsionSite" ) { p->torsionSite = stoi(param_val); }
    else if (input_param == "torsionSin2" ) { p->torsionSin2 = stoi(param_val); }
    else if (input_param == "torsionSin2V0" ) { p->torsionSin2V0 = stof(param_val); }
    else if (input_param == "torsionSin2V1" ) { p->torsionSin2V1 = stof(param_val); }
    else if (input_param == "torsionSin2omega" ) { p->torsionSin2omega = stof(param_val); }
    else if (input_param == "torsionSin2phi" ) { p->torsionSin2phi = stof(param_val); }
    else {  }
    getline (bash_in,line);
  }

  // close input file
  bash_in.close();

  // Error checking
  if ((p->VBPopFlag && p->CBPopFlag) || (p->VBPopFlag && p->QDPopFlag) || (p->CBPopFlag && p->QDPopFlag)) {
    std::cerr << "\nWARNING: population starting in multiple locations.\n";
  }

  if (p->CBPopFlag == POP_CONSTANT) {
    if (p->Nk_final > p->Nk || p->Nk_final < 1) {
      fprintf(stderr, "ERROR [Inputs]: Nk_final greater than Nk or less than 1.\n");
      exit(-1);
    }
    if (p->Nk_final < p->Nk_first) {
      fprintf(stderr, "ERROR [Inputs]: Nk_final is less than Nk_first.\n");
      exit(-1);
    }
  }

  if (p->QDPopFlag == POP_CONSTANT) {
    if (p->Nc_final > p->Nc || p->Nc_final < 1) {
      fprintf(stderr, "ERROR [Inputs]: Nc_final greater than Nc or less than 1.\n");
      exit(-1);
    }
    if (p->Nc_final < p->Nc_first) {
      fprintf(stderr, "ERROR [Inputs]: Nc_final is less than Nc_first.\n");
      exit(-1);
    }
  }

  if (p->Nl < 0) {
    fprintf(stderr, "ERROR [Inputs]: Nl less than 0.\n");
    exit(-1);
  }

  if (p->random_seed < -1) {
    std::cerr << "\nERROR: random_phase must be -1 or greater.\n";
    exit(-1);
  }

  if ((p->Nk < 2) && (p->rta)) {
    std::cerr << "\nERROR: when using RTA it is better to have many states in the conduction band." << std::endl;
  }

  // Decide which output files to make /////////////////////////////////////////

#ifdef DEBUG
  std::cout << "Assigning outputs as specified in " << p->inputFile << "\n";
#endif
  assignOutputs(p->inputFile.c_str(), p->outs, p);

  // error checking on various parameters //////////////////////////////////////

  // check torsion parameters, set up torsion spline ///////////////////////////

  if (p->torsion) {
#ifdef DEBUG
    std::cout << "Torsion is on." << std::endl;
#endif

    // error checking
    if (p->torsionSite > p->Nb) {
      std::cerr << "ERROR: torsion site (" << p->torsionSite
        << ") is larger than number of bridge sites (" << p->Nb << ")." << std::endl;
      exit(-1);
    }
    else if (p->torsionSite < 0) {
      std::cerr << "ERROR: torsion site is less than zero." << std::endl;
      exit(-1);
    }

    if (!p->torsionSin2) {
      if (!fileExists(p->torsionFile)) {
      std::cerr << "ERROR: torsion file " << p->torsionFile << " does not exist." << std::endl;
      exit(-1);
      }

      // create spline
      p->torsionV.readFile(p->torsionFile.c_str());
      if (p->torsionV.getFirstX() != 0.0) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should start at 0.0." << std::endl;
        exit(-1);
      }

      if (p->torsionV.getLastX() < p->tout) {
        std::cerr << "ERROR: time in " << p->torsionFile << " should be >= tout." << std::endl;
        exit(-1);
      }
    }
  }

  return;
}

// checks if a file exists (can be opened)
bool fileExists(std::string fileName) {
  std::ifstream ifile(fileName.c_str());
  return ifile.good();
}
