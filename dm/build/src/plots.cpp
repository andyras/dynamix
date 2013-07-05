#include "plots.hpp"

void makePlots(std::map<std::string, bool> &outs, struct PARAMETERS * p) {
 // probabilities in k states
 if (outs["kprobs.plt"] && (p->Nk > 1)) {
  plotKProbs("kprobs.plt", p);
 }
 // probabilities in c states
 if (outs["cprobs.plt"] && (p->Nc > 1)) {
  plotCProbs(p);
 }
 return;
}

void plotKProbs(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
 std::cout << "\nMaking kprobs.plt" << std::endl;
#endif
 std::ofstream output(fileName);

 output << "#!/usr/bin/env gnuplot" << std::endl
        << std::endl
        << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << std::endl
        << std::endl
        << "set output '/dev/null'" << std::endl
        << "plot '<cut -d \" \" -f 2- ./outs/kprobs.out' u ($2*"
	<< p->tout << "/" << p->numOutputSteps << "):($1*(" 
	<< p->kBandTop << "-" << p->kBandEdge << ")/("
	<< p->Nk << "-1)):3 matrix with image" << std::endl
        << std::endl
        << "set output 'figures/kprobs.pdf'" << std::endl
        << "set title 'Electron probability density in bulk conduction band'" << std::endl
        << "unset key " << std::endl
        //<< "unset colorbox" << std::endl
        << "set border 0" << std::endl
        << "set tics scale 0" << std::endl
        << "set ylabel 'Energy above band edge (a.u.)'" << std::endl
        << "set xlabel 'Time (a.u.)'" << std::endl
        << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << std::endl
        << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << std::endl
        << std::endl
        << "set palette defined(\\" << std::endl
        << "0.0    0.8667  0.8667  0.8667,\\" << std::endl
        << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl
        << "0.125  0.9255  0.8275  0.7725,\\" << std::endl
        << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl
        << "0.25   0.9608  0.7686  0.6784,\\" << std::endl
        << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl
        << "0.375  0.9686  0.6941  0.5804,\\" << std::endl
        << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl
        << "0.5    0.9569  0.6039  0.4824,\\" << std::endl
        << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl
        << "0.625  0.9255  0.4980  0.3882,\\" << std::endl
        << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl
        << "0.75   0.8706  0.3765  0.3020,\\" << std::endl
        << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl
        << "0.875  0.7961  0.2431  0.2196,\\" << std::endl
        << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl
        << "1      0.7059  0.0157  0.1490\\" << std::endl
        << ")" << std::endl
        << std::endl
        << "replot" << std::endl;
 return;
}

/* Makes a gnuplot file to plot the QD populations over time */
void plotCProbs(struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
 std::cout << "\nMaking cprobs.plt" << std::endl;
#endif
 std::ofstream output("cprobs.plt");
 output << "#!/usr/bin/env gnuplot" << std::endl
        << std::endl
        << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << std::endl
        << std::endl
        << "set output '/dev/null'" << std::endl
        << "plot '<cut -d \" \" -f 2- ./outs/cprobs.out' u ($2*"
	<< p->tout << "/" << p->numOutputSteps << "):($1):3 matrix with image" << std::endl
        << std::endl
        << "set output 'figures/cprobs.pdf'" << std::endl
        << "set title 'Electron probability density in QD'" << std::endl
        << "unset key " << std::endl
        //<< "unset colorbox" << std::endl
        << "set border 0" << std::endl
        << "set tics scale 0" << std::endl
        << "set ylabel 'State (index above band edge)'" << std::endl
        << "set xlabel 'Time (a.u.)'" << std::endl
        << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << std::endl
        << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << std::endl
        << std::endl
        << "set palette defined(\\" << std::endl
        << "0.0    0.8667  0.8667  0.8667,\\" << std::endl
        << "0.0625 0.8980  0.8471  0.8196,\\" << std::endl
        << "0.125  0.9255  0.8275  0.7725,\\" << std::endl
        << "0.1875 0.9451  0.8000  0.7255,\\" << std::endl
        << "0.25   0.9608  0.7686  0.6784,\\" << std::endl
        << "0.3125 0.9686  0.7333  0.6275,\\" << std::endl
        << "0.375  0.9686  0.6941  0.5804,\\" << std::endl
        << "0.4375 0.9686  0.6510  0.5294,\\" << std::endl
        << "0.5    0.9569  0.6039  0.4824,\\" << std::endl
        << "0.5625 0.9451  0.5529  0.4353,\\" << std::endl
        << "0.625  0.9255  0.4980  0.3882,\\" << std::endl
        << "0.6875 0.8980  0.4392  0.3451,\\" << std::endl
        << "0.75   0.8706  0.3765  0.3020,\\" << std::endl
        << "0.8125 0.8353  0.3137  0.2588,\\" << std::endl
        << "0.875  0.7961  0.2431  0.2196,\\" << std::endl
        << "0.9375 0.7529  0.1569  0.1843,\\" << std::endl
        << "1      0.7059  0.0157  0.1490\\" << std::endl
        << ")" << std::endl
        << std::endl
        << "replot" << std::endl;
 /*
 output << "#!/usr/bin/env gnuplot\n\n"
 << "reset\n"
 << "set terminal pdfcairo enhanced size 4in,3in font 'Arial-Bold,14'\n"
 << "set output '/dev/null'\n"
 << "!transpose -o _transpose ../outs/cprobs.out\n"
 << "plot '../outs/cprobs_transpose.out' every :::1 u ($1*" << p->tout << "/" << p->numOutputSteps << "):(-$2):3 matrix with image\n"
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
 */

 return;
}
