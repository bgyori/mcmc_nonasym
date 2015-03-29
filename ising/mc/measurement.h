#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <statistics.h>
#include <vector>

struct measurement {
	measurement() {saveFx=0;}
	measurement(int savef) {saveFx = savef;}
		void accumulate_ratio(int a,double T) {
		ratio.accumulate((double)a,T);
		}
	void accumulate_fx(double f,double T) {
		fx.accumulate(f,T);
		if(saveFx) fxAll.push_back(f);
		}

	void getFxAll(std::vector<double>& vIn){
		vIn.swap(fxAll);
		}
	
	template<class outstream>
	void print_ratio(outstream& out) const { ratio.print(out); }
	template<class outstream>
	void print_fx(outstream& out) const { fx.print(out); }

	statistics ratio,fx;
	std::vector<double> fxAll;
	
	private:
		int saveFx;
	};

#endif

