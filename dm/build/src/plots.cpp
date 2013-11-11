#include "plots.hpp"

//#define DEBUG_PLOT

void makePlots(std::map<const std::string, bool> &outs, struct PARAMETERS * p) {
  // populations in subsystems
  try {
    if (outs.at("populations.plt")) {
      plotPopulations("populations.plt", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // probabilities in k states
    if (outs.at("kprobs.plt") && (p->Nk > 1)) {
      plotKProbs("kprobs.plt", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // probabilities in c states
    if (outs.at("cprobs.plt") && (p->Nc > 1)) {
      plotCProbs(p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // density matrix in time
    if (outs.at("dmt_z.plt")) {
      plotDMt_z("dmt_z.plt", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  try {
    // populations in k states as a movie
    if (outs.at("kprobs_movie.plt") && (p->Nk > 1)) {
      plotKProbsMovie("kprobs_movie.plt", p);
    }
  }
  catch (const std::out_of_range& oor) {
#ifdef DEBUG
    std::cerr << "Out of Range error: " << oor.what() << std::endl;
#endif
  }

  return;
}

/* Plots the populations in different subsystems over time,
 * separately and on the same plot. */
void plotPopulations(char * fileName, struct PARAMETERS * p) {
  std::ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << std::endl;
  o << std::endl;
  o << "reset" << std::endl;
  o << std::endl;
  o << "set style data lines" << std::endl;
  o << "set style line 1 lt 1 lc rgb 'red'" << std::endl;
  o << "set style line 2 lt 1 lc rgb 'dark-green'" << std::endl;
  o << "set style line 3 lt 1 lc rgb 'blue'" << std::endl;
  o << "set style increment user" << std::endl;
  o << std::endl;
  o << "set terminal pdfcairo enhanced dashed size 5,5 lw 2 font 'Arial-Bold,12'" << std::endl;
  o << "set output 'figures/populations.pdf'" << std::endl;
  o << std::endl;
  o << "stats './outs/tkprob.out' u ($1/41.3414):2 nooutput" << std::endl;
  o << "set xrange [STATS_min_x:STATS_max_x]" << std::endl;
  o << std::endl;
  o << "set rmargin 3" << std::endl;
  o << "set tmargin 0" << std::endl;
  o << "set bmargin 3" << std::endl;
  o << "set tics scale 0 nomirror" << std::endl;
  o << "unset xtics" << std::endl;
  o << "set format y '%.2e'" << std::endl;
  o << std::endl;
  o << "set multiplot" << std::endl;
  if (p->bridge_on) {
    o << "set size 1,0.31" << std::endl;
  }
  else {
    o << "set size 1,0.47" << std::endl;
  }
  o << std::endl;
  if (p->bridge_on) {
    o << "set origin 0,0.62" << std::endl;
  }
  else {
    o << "set origin 0,0.48" << std::endl;
  }
  o << "set ylabel 'Population (a.u.)'" << std::endl;
  o << "set title 'Bulk Population vs. Time'" << std::endl;
  o << "plot './outs/tkprob.out' u ($1/41.3414):2 lt 1 notitle" << std::endl;
  o << std::endl;
  if (p->bridge_on) {
    o << "set origin 0,0.32" << std::endl;
    o << "set title 'Bridge Population vs. Time'" << std::endl;
    o << "plot './outs/bprobs.out' u ($1/41.3414):2 lt 3 notitle";
    if (p->Nb > 1) {
      for (int ii = 0; ii < (p->Nb - 1); ii++) {
	o << ", \\" << std::endl << "     '' u ($1/41.3414):" << (ii+3) << " lt " << (ii+4) << " notitle";
      }
    }
    o << std::endl << std::endl;
  }
  o << "set origin 0,0.01" << std::endl;
  o << "set xlabel 'Time (fs)'" << std::endl;
  o << "set xtics scale 0" << std::endl;
  o << "set title 'QD Population vs. Time'" << std::endl;
  o << "plot './outs/tcprob.out' u ($1/41.3414):2 lt 2 notitle" << std::endl;
  o << std::endl;
  o << "unset multiplot" << std::endl;
  o << std::endl;
  o << "set terminal pdfcairo enhanced size 5,3 lw 2 font 'Arial-Bold,12' dashed rounded" << std::endl;
  o << "set output 'figures/populationsTogether.pdf'" << std::endl;
  o << std::endl;
  o << "set key tmargin right" << std::endl;
  o << "set size 1,1" << std::endl;
  o << "set tmargin 4" << std::endl;
  o << "set ytics 1 format '%.0f'" << std::endl;
  o << "unset style line" << std::endl;
  o << "set style line 1 lc rgb 'red'" << std::endl;
  o << "set style line 2 lc rgb 'dark-green'" << std::endl;
  o << "set style line 3 lc rgb 'blue'" << std::endl;
  o << std::endl;
  o << "set title 'Subsystem Populations vs. Time' offset 0,0.6" << std::endl;
  o << "plot './outs/tkprob.out' u ($1/41.3414):2 lt 1 lw 2 title 'Bulk', \\" << std::endl;
  if (p->bridge_on) {
    o<< "'./outs/tbprob.out' u ($1/41.3414):2 lt 3 lw 2 title 'Bridge', \\" << std::endl;
  }
  o<< "'./outs/tcprob.out' u ($1/41.3414):2 lt 2 lw 2 title 'QD'" << std::endl;

  return;
}

void plotKProbs(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
  std::cout << "\nMaking kprobs.plt" << std::endl;
#endif
  std::ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << std::endl;
  o << std::endl;
  o << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << std::endl;
  o << std::endl;
  o << "set output '/dev/null'" << std::endl;
  o << "plot '<cut -d \" \" -f 2- ./outs/kprobs.out' u ($2*"
    << p->tout << "/" << p->numOutputSteps << "):($1*(" 
    << p->kBandTop << "-" << p->kBandEdge << ")/(";
  o << p->Nk << "-1)):3 matrix with image" << std::endl;
  o << std::endl;
  o << "set output 'figures/kprobs.pdf'" << std::endl;
  o << "set title 'Electron probability density in bulk conduction band'" << std::endl;
  o << "unset key " << std::endl;
  o << "set border 0" << std::endl;
  o << "set tics scale 0" << std::endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << std::endl;
  o << "set xlabel 'Time (a.u.)'" << std::endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << std::endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << std::endl;
  o << std::endl;
  o << "set palette defined(\\" << std::endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << std::endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << std::endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << std::endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << std::endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << std::endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << std::endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << std::endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << std::endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl;
  o << "1      0.7059  0.0157  0.1490\\" << std::endl;
  o << ")" << std::endl;
  o << std::endl;
  o << "replot" << std::endl;

  return;
}

/* Makes a gnuplot file to plot the QD populations over time */
void plotCProbs(struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
  std::cout << "\nMaking cprobs.plt" << std::endl;
#endif
  std::ofstream o("cprobs.plt");

  o << "#!/usr/bin/env gnuplot" << std::endl;
  o << std::endl;
  o << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << std::endl;
  o << std::endl;
  o << "set output '/dev/null'" << std::endl;
  o << "plot '<cut -d \" \" -f 2- ./outs/cprobs.out' u ($2*"
    << p->tout << "/" << p->numOutputSteps << "):($1):3 matrix with image" << std::endl;
  o << std::endl;
  o << "set output 'figures/cprobs.pdf'" << std::endl;
  o << "set title 'Electron probability density in QD'" << std::endl;
  o << "unset key " << std::endl;
  o << "set border 0" << std::endl;
  o << "set tics scale 0" << std::endl;
  o << "set ylabel 'State (index above band edge)'" << std::endl;
  o << "set xlabel 'Time (a.u.)'" << std::endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << std::endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << std::endl;
  o << std::endl;
  o << "set palette defined(\\" << std::endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << std::endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << std::endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << std::endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << std::endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << std::endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << std::endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << std::endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << std::endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl;
  o << "1      0.7059  0.0157  0.1490\\" << std::endl;
  o << ")" << std::endl;
  o << std::endl;
  o << "replot" << std::endl;

  return;
}

/* Plots the density matrix in time as a series of .png files */
void plotDMt_z(char * fileName, struct PARAMETERS * p) {
  std::ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << std::endl;
  o << std::endl;
  o << "![ -d img ] && rm -rf img" << std::endl;
  o << "!mkdir -p img" << std::endl;
  o << std::endl;
  o << "set palette defined(\\" << std::endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << std::endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << std::endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << std::endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << std::endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << std::endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << std::endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << std::endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << std::endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl;
  o << "1      0.7059  0.0157  0.1490\\" << std::endl;
  o << ")" << std::endl;
  o << std::endl;
  o << "set terminal pngcairo font 'Arial-Bold,12'" << std::endl;
  o << "set output '/dev/null'" << std::endl;
  o << "plot 'outs/dmt_z.out' index 0 matrix with image" << std::endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << std::endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << std::endl;
  o << std::endl;
  o << "set tics out scale 0.5 nomirror" << std::endl;
  o << "set size square" << std::endl;
  o << std::endl;
  o << "do for [ii=0:" << p->numOutputSteps << "] {" << std::endl;
  o << "set output sprintf(\"img/dm%.5d.png\", ii)" << std::endl;
  o << "set title sprintf(\"Density Matrix at t = %f fs\", 1.0*ii*"
    << p->tout << "/" << p->numOutputSteps << ")" << std::endl;
  o << "plot 'outs/dmt_z.out' index ii matrix with image" << std::endl;
  o << "}" << std::endl;
  o << std::endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o /dev/null \"mf://img/*.png\"" << std::endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o figures/dmt_z.avi \"mf://img/*.png\"" << std::endl;
  o << "!rm -f divx2pass.log" << std::endl;
  o << "!rm -rf img" << std::endl;
  return;
}

/* Plots the populations in the bulk states over time as a movie */
void plotKProbsMovie(char * fileName, struct PARAMETERS * p) {
  std::ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << std::endl;
  o << std::endl;
  o << "![ -d img ] && rm -rf img" << std::endl;
  o << "!mkdir -p img" << std::endl;
  o << "!transpose -o _t outs/kprobs.out" << std::endl;
  o << std::endl;
  o << "set palette defined(\\" << std::endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << std::endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << std::endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << std::endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << std::endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << std::endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << std::endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << std::endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << std::endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl;
  o << "1      0.7059  0.0157  0.1490\\" << std::endl;
  o << ")" << std::endl;
  o << std::endl;
  o << "set terminal pngcairo font 'Arial-Bold,12'" << std::endl;
  o << "set tics scale 0" << std::endl;
  o << "unset key" << std::endl;
  o << "do for [ii=0:" << p->numOutputSteps << "] {" << std::endl;
  o << "idx=ii+1" << std::endl;
  o << "set output sprintf(\"img/%.5d.png\", ii)" << std::endl;
  o << "set multiplot layout 2,1" << std::endl;
  o << std::endl;
  o << "set xrange [0:" << p->kBandWidth << "]" << std::endl;
  o << "set yrange [0:*]" << std::endl;
  o << "set xlabel 'Energy above band edge (a.u.)'" << std::endl;
  o << "set ylabel 'Population in state'" << std::endl;
  o << "set title sprintf(\"Bulk Populations at t = %.6f a.u.\", ii*"
    << p->tout << "/" << p->numOutputSteps << ")" << std::endl;
  o << "plot 'outs/kprobs_t.out' u ($0*" << p->kBandWidth << "/" << (p->Nk - 1)
    << "):idx every ::1 with filledcurves x1" << std::endl;
  o << std::endl;
  o << "set xrange [0:" << p->tout << "]" << std::endl;
  o << "set yrange [0.00:" << p->kBandWidth << "]" << std::endl;
  o << "set xlabel 'Time (a.u.)'" << std::endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << std::endl;
  o << "set title 'State Populations vs. Energy and Time'" << std::endl;
  o << "plot '<cut -d \" \" -f 2- outs/kprobs.out' u ($2*" << p->tout << "/"
    << p->numOutputSteps << "):($1*" << p->kBandWidth << "/" << (p->Nk - 1) << "):3 matrix with image" << std::endl;
  o << "unset multiplot" << std::endl;
  o << "}" << std::endl;
  o << std::endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o /dev/null \"mf://img/*.png\"" << std::endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o figures/kprobs_movie.avi \"mf://img/*.png\"" << std::endl;
  o << "!rm -f divx2pass.log" << std::endl;
  o << "!rm -rf img" << std::endl;
  return;
}
