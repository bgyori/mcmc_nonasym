#ifndef STATISTICS_H
#define STATISTICS_H

#include <map>

struct _statData {
  double sum1;
  double sum2;
  double num;
  double mean;
  double stddev;
  double stde;
};

struct statistics {

  typedef _statData value_type;

  statistics(){}
  
  void accumulate(double a,double T) {
    if(data.find(T)==data.end()) {
      value_type v;
      v.sum1 = a; v.sum2 = a*a; v.num = 1; v.mean = a;
      v.stddev = 0;
      data.insert(std::make_pair(T,v));
    }
    else { 
      data[T].sum1  += a; 
      data[T].sum2  += a*a; 
      data[T].num   += 1;
      double m1 = data[T].sum1/data[T].num;
      double m2 = data[T].sum2/data[T].num;
      data[T].mean   = m1;
      data[T].stddev = sqrt(m2-m1*m1);
      data[T].stde = data[T].stddev / sqrt(data[T].num-1);
    }
  }
  
  template<class outstream>
  void print(outstream& out) const {

    std::map<double,value_type  >::const_iterator it=data.begin();
    
    for(;it!=data.end();++it) {
      out<<it->first<<" "<<it->second.mean<<" "<<it->second.stddev
	 <<" "<<it->second.stde<<std::endl;
    }
  }
  
  template<class outstream>
  void printMean(outstream& out) const {
    std::map<double,value_type>::const_iterator it=data.begin();
    for(;it!=data.end();++it) {
      out<<it->second.mean<<std::endl;
      }
    }
  
  
  double getMean() const {
	return ((data.begin())->second).mean;
	}
  
  double getVar() const {
	double s1 = (data.begin()->second).sum1;
	double s2 = (data.begin()->second).sum2;
	double n = (data.begin()->second).num;
	return s2-s1*(s1/n);
  }

  private:
    std::map<double,value_type>     data;
};

#endif

