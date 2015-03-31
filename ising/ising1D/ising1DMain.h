#include <isingLattice1D.h>
#include "hamiltonian.h"

#include "glauber_rate.h"
#include "measurement.h"
#include "rng_mersenne.h"
#include "monte_carlo.h"

#include <iostream>
#include <fstream>

#ifdef METROPOLIS
	#include "single_site_update.h"
	#include "metropolis_rate.h"
	typedef single_site_update<isingLattice1D,hamiltonian,metropolis_rate> updater_type;
#else
	#include "glauber_rate.h"
	#ifdef SWEEP
		#include "sweep_update.h"
		typedef sweep_update<isingLattice1D,hamiltonian,glauber_rate> updater_type;
	#else
		#include "single_site_update.h"
		typedef single_site_update<isingLattice1D,hamiltonian,glauber_rate> updater_type;
	#endif
#endif
	
typedef monte_carlo<updater_type,isingLattice1D,measurement>            mc_type;

void magIsing1D(int nRuns, int latticeSize, int nSteps, int nRelaxSteps, int nStepInterval, 
		double T, double H, int seed, double* magMean, double* magVar, std::vector<double>& mags){
  
  hamiltonian                   energy(H,1.0);
  #ifdef METROPOLIS
		metropolis_rate				rate;
	#else
		glauber_rate				rate;
	#endif
  
  rng_mersenne ran(seed);

  
  if(!magMean){
	isingLattice1D 	lattice(latticeSize,ran);
	updater_type        updater(lattice,energy,rate);
	measurement         data(1);
	mc_type             mc(updater,lattice,data);

	mc(ran,nRelaxSteps,T);
	mc(ran,nSteps,T,nStepInterval);
	data.getFxAll(mags);
	}
  else{
	  for (int i=0;i<nRuns;i++){
		isingLattice1D 	lattice(latticeSize,ran);
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
