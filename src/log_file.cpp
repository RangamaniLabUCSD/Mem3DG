#include "ddgsolver/force.h"
#include "ddgsolver/integrator.h"

#include <iostream>
#include <fstream>
using namespace std;

namespace ddgsolver {
	void integrator::getLogFiles() {
		ofstream myfile("output-file/output.txt");
		if (myfile.is_open())
		{

			myfile << "Physical parameters used: \n";
			myfile << "\n";
			myfile << "Kb:     " << f.Kb << "\n"
				<< "H0:     " << f.H0 << "\n"
				<< "Kse:    " << f.Kse << "\n"
				<< "Ksl:    " << f.Ksl << "\n"
				<< "Ksg:    " << f.Ksg << "\n"
				<< "Kv:     " << f.Kv << "\n"
				<< "gamma:  " << f.gamma << "\n"
				<< "Vt:     " << f.Vt << "\n"
				<< "kt:     " << f.kt << "\n"
				<< "sigma:  " << f.sigma << "\n"
				<< "ptInd:  " << f.ptInd << "\n"
				<< "extF:   " << f.extF << "\n"
				<< "conc:   " << f.conc << "\n";
			myfile << "\n";
			myfile << "Integration parameters used: \n";
			myfile << "\n";
			myfile << "dt:    " << this -> dt << "\n"
						 << "T:     " << this -> total_time << "\n"
						 << "eps:		" << this -> tolerance << "\n"
						 << "tSave: " << this -> tSave << "\n";

			myfile.close();
		}
		else cout << "Unable to open file";
	}
}
