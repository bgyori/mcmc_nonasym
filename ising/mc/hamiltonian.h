#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <iostream>

struct hamiltonian {
  hamiltonian() {H = 0;J = 1;}
  hamiltonian(double h,double j): H(h), J(j) {}
  
  template<class model_type>
  double operator()(model_type& ising) const {
 
    typedef typename model_type::iterator iterator;
    double ret = 0;
    int N = ising.size();
    for(iterator it=ising.begin();it!=ising.end();++it) {
      ret += -J*(*it)*ising.effective_field(it);
	  ret += -H*(*it);
	}
    return ret/(1.0*N);
  }

  template<class model_type>
  double operator()(model_type& ising,
                    typename model_type::iterator& it,
                    typename model_type::value_type& new_value) const {
    double ret;
    ret = -J*(new_value-(*it)) * ising.effective_field(it);
    ret += -H*(new_value-(*it));
	return ret;
  }
  
private:
	double H,J;
};

#endif

