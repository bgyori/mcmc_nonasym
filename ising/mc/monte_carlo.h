#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H
#include <iostream>

template<class updater_type,
         class model_type,
	 class data_type>
struct monte_carlo {

  monte_carlo(updater_type& u,
              model_type& m,
              data_type& d) 
    : updater(u),model(m),data(d) { }

  template<class random_type>
  void operator()(random_type& ran, unsigned int mcs,
                  double T,unsigned int interval) {

	unsigned int acc;
	for(unsigned int i=0;i<mcs;i++){
		acc = updater(ran,T);
		data.accumulate_ratio(acc,T);
		if(i%interval == 0){
			#ifdef MAGSIGN
				data.accumulate_fx(model.magSign(),T);
			#else
				data.accumulate_fx(model.magnetization(),T);
			#endif
			}
		}
  }

template<class random_type>
  void operator()(random_type& ran, unsigned int mcs,
                  double T) {

	for(unsigned int i=0;i<mcs;i++){
		updater(ran,T);
		}	
  }
 
  private:
    updater_type&    updater;
    model_type&      model;
    data_type&       data;
};

#endif

