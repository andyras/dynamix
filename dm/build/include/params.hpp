#pragma once

#include <vector>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>
#include <nvector/nvector_serial.h>
#include <mkl.h>

#include "constants.hpp"
#include "spline.hpp"
#include "numerical.hpp"

#ifdef __BOOST_SERIALIZE__
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#endif

class Params {
public:
  int nproc = 0;                        // number of processors
  bool wavefunction = 1;                // 1 => wavefunction; 0 => density matrix
  bool coherent = 1;                    // 1 => \rho = \ket{\psi}\bra{\psi}; 0 => \rho_{ii} only
  bool justPlots = false;               // just make plots, no propagation or other output
  bool timedepH = true;                 // if H is TD, use CVODE, else diag H and propogate
  bool analytical = false;              // turn on analytical propagation
  bool kinetic = false;                 // kinetic relaxation model
  bool kineticQD = false;               // kinetic relaxation model in QD
  bool dynamicMu = false;               // for kinetic model, calculate Fermi level dynamically
  bool dynamicMuQD = false;             // for kinetic model, calculate Fermi level dynamically on QD
  bool rta = true;                      // turn on relaxation time approximation (RTA)
  bool rtaQD = true;                    // turn on relaxation time approximation (RTA)
  bool dephasing = false;               // turn on dephasing
  bool progressFile = false;            // create a file to show progress of the run
  bool progressStdout = false;          // show progress in stdout
  realtype abstol = 1e-10;              // absolute tolerance (for SUNDIALS)
  realtype reltol = 1e-10;              // relative tolerance (for SUNDIALS)
  realtype tout = 10000;                // final time reached by solver in atomic units
  int numsteps = 10000;                 // number of time steps
  int numOutputSteps = 1000;            // number of timesteps
  int NEQ = 1;                          // total number of states
  int NEQ2 = 1;                         // total number of coefficients in density matrix
  realtype kBandEdge = 0.0;             // lower edge of bulk conduction band
  realtype kBandTop = 0.01;             // upper edge of bulk conduction band
  realtype lBandTop = -0.01;            // upper edge of bulk valence band
  realtype bulk_gap = 0.001;            // bulk band gap
  double valenceBand = 0.01;            // valence band width
  double bulkGaussSigma = 0.001;        // width of initial Gaussian in bulk
  double bulkGaussMu = 0.01;            // position of initial Gaussian above band edge
  double me = 1.0;                      // effective mass of electron
  double mh = 1.0;                      // effective mass of hole
  double me_c = 1.0;                    // effective mass of electron on QD
  double mh_c = 1.0;                    // effective mass of hole on QD
  double X2 = 1512.2873345935727;       // "Bohr radius" of material, inverse spacing in k-space
  realtype temperature = 3e2;           // temperature of the system
  double EF = 0.0;                      // Fermi level in bulk
  realtype gamma1 = 1e-3;               // \gamma_1 in RTA (relaxation rate)
  realtype gamma1_c = 1e-3;             // \gamma_1 in RTA (relaxation rate) on QD
  realtype gamma2 = 1e-3;               // \gamma_2 in RTA (dephasing rate)
  double muLK = 1.0;                    // transition dipole moment from l to k (energy a.u.)
  double pumpFWHM = 1000;               // FWHM of pump pulse (time a.u.)
  double pumpPeak = 2000;               // time of peak of pump pulse (a.u.)
  double pumpFreq = 0.01;               // frequency of pump pulse (energy a.u.)
  double pumpAmpl = 1.0;                // intensity of pump pulse (electric field a.u.)
  double pumpPhase = 0.0;               // pump pulse phase (in units of radians)
  int CBPopFlag = 0;                    // flag for starting condition in conduction band
  int VBPopFlag = 0;                    // flag for starting condition in valence band
  int QDPopFlag = 0;                    // flag for starting condition in QD
  int Nk_first = 1;     // first k state initially populated
  int Nk_final = 1;     // final k state initially populated
  int Nc_first = 1;     // first c state initially populated
  int Nc_final = 1;     // final c state initially populated
  bool laser_on = false;
  bool parabolicCoupling = false;
  bool scale_bubr = false;
  bool scale_brqd = false;
  bool scale_buqd = false;
  bool scale_laser = false;
  bool bridge_on = false;
  bool random_phase = false;
  int random_seed = 0;
  bool torsion = false;
  std::string torsionFile;
  int torsionSite = 0;
  bool torsionSin2 = false;             // flag to turn on sin^2 torsional coupling
  bool torsionGaussianPulse = false;    // flag for Gaussian "half-cycle" torsion coupling
  bool torsionCos2Pulse = false;        // flag for cos^2 "half-cycle" torsion coupling
  double torsionCouplingV0 = 0.001;     // V(t) = V0 + V1*sin^2(omega*t + phi)
  double torsionCouplingV1 = 0.001;
  double torsionCouplingOmega = 0.001;
  double torsionCouplingPhi = 0.0;

  int Nk = 1;                           // number of each type of state
  int Nc = 0;                           // number of each type of state
  int Nb = 0;                           // number of each type of state
  int Nl = 0;                           // number of each type of state
  int Ik = 0;                           // index starters for each type of state
  int Ic = 0;                           // index starters for each type of state
  int Ib = 0;                           // index starters for each type of state
  int Il = 0;                           // index starters for each type of state

  realtype kBandWidth;

  std::vector<realtype> energies;
  std::vector<realtype> k_energies;
  std::vector<realtype> c_energies;
  std::vector<realtype> b_energies;
  std::vector<realtype> l_energies;
  std::vector<realtype> Vbridge;
  std::vector<realtype> Vnobridge;
  std::vector< std::vector<realtype> > V;
  std::vector<realtype> H;
  std::vector<realtype> H_sp;
  std::vector<int> H_cols;
  std::vector<int> H_rowind;
  std::vector<realtype> times;
  std::vector<realtype> startWfn;
  std::vector<realtype> startDM;
  std::vector<realtype> wfnt;
  std::vector<realtype> dmt;

  // map of output file names to bool
  std::map<const std::string, bool> outs;

  // output directory
  std::string outputDir = "outs/";

  // input file names
  std::string inputFile = "ins/parameters.in";
  std::string cEnergiesInput = "ins/c_energies.in";
  std::string bEnergiesInput = "ins/b_energies.in";
  std::string VNoBridgeInput = "ins/Vnobridge.in";
  std::string VBridgeInput = "ins/Vbridge.in";

  Spline torsionV;

  double lastTime;                      // value of most recently calculated timepoint

  double lastMu;                        // value of Fermi level at last time point
  double lastMuQD;                      // value of Fermi level in QD at last time point

  // methods ///////////////////////////////////////////////////////////////////
  void buildHamiltonian();

  void buildCoupling();

  double getTorsionCoupling(double t);

#ifdef __BOOST_SERIALIZE__
private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & nproc;
    ar & wavefunction;
    ar & coherent;
    ar & justPlots;
    ar & timedepH;
    ar & analytical;
    ar & kinetic;
    ar & kineticQD;
    ar & dynamicMu;
    ar & dynamicMuQD;
    ar & rta;
    ar & rtaQD;
    ar & dephasing;
    ar & progressFile;
    ar & progressStdout;
    ar & abstol;
    ar & reltol;
    ar & tout;
    ar & numsteps;
    ar & numOutputSteps;
    ar & NEQ;
    ar & NEQ2;
    ar & kBandEdge;
    ar & kBandTop;
    ar & lBandTop;
    ar & bulk_gap;
    ar & valenceBand;
    ar & bulkGaussSigma;
    ar & bulkGaussMu;
    ar & me;
    ar & mh;
    ar & me_c;
    ar & mh_c;
    ar & X2;
    ar & temperature;
    ar & EF;
    ar & gamma1;
    ar & gamma1_c;
    ar & gamma2;
    ar & muLK;
    ar & pumpFWHM;
    ar & pumpPeak;
    ar & pumpFreq;
    ar & pumpAmpl;
    ar & pumpPhase;
    ar & CBPopFlag;
    ar & VBPopFlag;
    ar & QDPopFlag;
    ar & Nk_first;
    ar & Nk_final;
    ar & Nc_first;
    ar & Nc_final;
    ar & laser_on;
    ar & parabolicCoupling;
    ar & scale_bubr;
    ar & scale_brqd;
    ar & scale_buqd;
    ar & scale_laser;
    ar & bridge_on;
    ar & random_phase;
    ar & random_seed;
    ar & torsion;
    ar & torsionFile;;
    ar & torsionSite;
    ar & torsionSin2;
    ar & torsionGaussianPulse;
    ar & torsionCos2Pulse;
    ar & torsionCouplingV0;
    ar & torsionCouplingV1;
    ar & torsionCouplingOmega;
    ar & torsionCouplingPhi;
    ar & Nk;
    ar & Nc;
    ar & Nb;
    ar & Nl;
    ar & Ik;
    ar & Ic;
    ar & Ib;
    ar & Il;

    ar & kBandWidth;

    ar & energies;
    ar & k_energies;
    ar & c_energies;
    ar & b_energies;
    ar & l_energies;
    ar & Vbridge;
    ar & Vnobridge;
    ar & V;
    ar & H;
    ar & H_sp;
    ar & H_cols;
    ar & H_rowind;
    ar & times;
    ar & startWfn;
    ar & startDM;
    ar & wfnt;
    ar & dmt;

    ar & outs;

    ar & outputDir;

    ar & inputFile;
    ar & cEnergiesInput;
    ar & bEnergiesInput;
    ar & VNoBridgeInput;
    ar & VBridgeInput;

    ar & torsionV;

    ar & lastTime;

    ar & lastMu;
    ar & lastMuQD;
  }
#endif
};