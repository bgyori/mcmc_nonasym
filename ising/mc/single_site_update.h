#ifndef SINGLE_SITE_UPDATE_H
#define SINGLE_SITE_UPDATE_H


#include <cmath>
#include <stdexcept>
#include <iostream>

template<class model_type,
         class hamiltonian_type,
         class acceptance_rate_type>
class single_site_update {

public:

  single_site_update(model_type& m, 
                     hamiltonian_type& e,
                     acceptance_rate_type& r)
    : model(m), energy(e), rate(r) {}

  template<class random_type>
  unsigned int operator()(random_type& ran,double T) {
	typename model_type::iterator it;
	typename model_type::value_type new_value;
	it = model.random_site(ran);			// Chosse random site
	new_value = -1*(*it);				// Flip spin as proposal
	double de = energy(model, it, new_value);	// Calculate energy change
	double rt = rate(de,T);				// Calculate probability of step
	double r = ran.randu();		// r=U[0,1]
	if(r < rt){					// Accept flip
		model.flip(it);
		return 1;
		}
	else {						// Reject flip
		return 0;
		}
  }

private:
  model_type&               model;
  hamiltonian_type&         energy;
  acceptance_rate_type&     rate;
};  


#endif 

