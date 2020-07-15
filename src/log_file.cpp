#include "ddgsolver/force.h"
#include "ddgsolver/integrator.h"

#include <fstream>
#include <iostream>
using namespace std;

namespace ddgsolver {
namespace integration {

void getLogFiles(Force &f, double dt, double total_time, double tolerance,
                 double tSave, std::string outputDir) {
  ofstream myfile(outputDir + "output.txt");
  if (myfile.is_open()) {

    myfile << "Physical parameters used: \n";
    myfile << "\n";
    myfile << "Kb:     " << f.P.Kb << "\n"
           << "H0:     " << f.P.H0 << "\n"
           << "Kse:    " << f.P.Kse << "\n"
           << "Ksl:    " << f.P.Ksl << "\n"
           << "Ksg:    " << f.P.Ksg << "\n"
           << "Kv:     " << f.P.Kv << "\n"
           << "gamma:  " << f.P.gamma << "\n"
           << "Vt:     " << f.P.Vt << "\n"
           << "kt:     " << f.P.kt << "\n"
           << "sigma:  " << f.P.sigma << "\n"
           << "ptInd:  " << f.P.ptInd << "\n"
           << "extF:   " << f.P.extF << "\n"
           << "conc:   " << f.P.conc << "\n";
    myfile << "\n";
    myfile << "Integration parameters used: \n";
    myfile << "\n";
    myfile << "dt:    " << dt << "\n"
           << "T:     " << total_time << "\n"
           << "eps:		" << tolerance << "\n"
           << "tSave: " << tSave << "\n";

    myfile.close();
  } else
    cout << "Unable to open file";
}

} // namespace integration
} // namespace ddgsolver
