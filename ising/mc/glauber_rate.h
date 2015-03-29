#ifndef GLAUBER_RATE_H
#define GLAUBER_RATE_H

#include <cmath>
#include <iostream>
struct glauber_rate {

  double operator()(double de, double T) {
	double e = 1/(1+exp(de/T));
	return e;
  }
};

#endif
