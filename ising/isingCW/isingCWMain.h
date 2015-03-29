#include "isingLatticeCW.h"
#include "hamiltonian.h"
#include "single_site_update.h"

#include "measurement.h"
#include "sml_shift_register.h"
#include "monte_carlo.h"

#include <iostream>
#include <fstream>
#include <vector>

#ifdef METROPOLIS
	#include "metropolis_rate.h"
	typedef single_site_update<isingCW,hamiltonian,metropolis_rate> updater_type;
#else
	#include "glauber_rate.h"
	typedef single_site_update<isingCW,hamiltonian,glauber_rate> updater_type;
#endif

typedef monte_carlo<updater_type,isingCW,measurement>            mc_type;

void magIsingCW(int nRuns, int latticeSize, int nSteps, int nRelaxSteps, int nStepInterval, 
					double T, double H, int seed, double* magMean, double* magVar, std::vector<double>& mags){


	hamiltonian					energy(H,(1.0/latticeSize));
	
	#ifdef METROPOLIS
		metropolis_rate				rate;
	#else
		glauber_rate				rate;
	#endif
	sml::shift_register_default	ran(seed);

	if(!magMean){
		isingCW 	lattice(latticeSize,ran);
		updater_type        updater(lattice,energy,rate);
		measurement         data(1);
		mc_type             mc(updater,lattice,data);
		
		mc(ran,nRelaxSteps,T);
		mc(ran,nSteps,T,nStepInterval);
		data.getFxAll(mags);
		}
	else {
		for (int i=0;i<nRuns;i++){
			isingCW 	lattice(latticeSize,ran);
			updater_type        updater(lattice,energy,rate);
			measurement         data;
			mc_type             mc(updater,lattice,data);

			mc(ran,nRelaxSteps,T);
			mc(ran,nSteps,T,nStepInterval);
			magMean[i] = data.fx.getMean();
			magVar[i] = data.fx.getVar();
		}
	}
	

}


