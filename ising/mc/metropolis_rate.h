#ifndef METROPOLIS_RATE_H
#define METROPOLIS_RATE_H

#include <cmath>
#include <iostream>
struct metropolis_rate {

  double operator()(double de, double T) {
	if(de>0)
		return exp(-de/T); 
	else
		return 1;
  }
};

#endif
