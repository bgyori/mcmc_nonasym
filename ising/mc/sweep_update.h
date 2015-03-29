#ifndef SWEEP_UPDATE_H
#define SWEEP_UPDATE_H


#include <cmath>
#include <stdexcept>
#include <iostream>

template<class model_type,
         class hamiltonian_type,
         class acceptance_rate_type>
class sweep_update {

public:

  sweep_update(model_type& m, 
               hamiltonian_type& e,
               acceptance_rate_type& r)
    : model(m), energy(e), rate(r) {}

  template<class random_type>
  unsigned int operator()(random_type& ran,double T) {
	typename model_type::iterator it;
	typename model_type::value_type new_value;
	double r, rt, de;
	
	for(it=model.begin();it!=model.end();++it){
		new_value = -1*(*it);				// Flip spin as proposal
		de = energy(model, it, new_value);	// Calculate energy change
		rt = rate(de,T);					// Calculate probability of step
		r = (ran()*1.0 / ran.max());		// r=U[0,1]
		if(r < rt){							// Accept flip
			model.flip(it);
			}
		}
	return 1;
  }

private:
  model_type&               model;
  hamiltonian_type&         energy;
  acceptance_rate_type&     rate;
};  


#endif 

