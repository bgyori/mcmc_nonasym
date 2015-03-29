#ifndef LATTICE1D_H
#define LATTICE1D_H


#include <numeric>
#include <iostream>
#include <vector>

template<class data_container_type>
struct _ising_iterator {

  typedef _ising_iterator<data_container_type>      this_type;
  typedef typename data_container_type::iterator    data_iterator;
  typedef typename data_container_type::value_type  value_type;

  _ising_iterator() { }
  _ising_iterator(data_iterator dit,unsigned int id) :
    ditr(dit), indx(id) { }

  void           operator++()                    { ditr++; indx++; }
  void           operator--()                    { ditr--; indx--; }
  value_type&    operator*()                     { return *ditr;  }
  unsigned int&  index()                         { return indx;   }
  bool           operator!=(const this_type& n)  { return !this->operator==(n); }
  bool           operator==(const this_type& n)  { 
    return ( (n.ditr==ditr) && (n.indx==indx) ); 
  }
  void           operator=(const this_type& n)   { ditr=n.ditr; indx=n.indx; }
  private:
    data_iterator   ditr;
    unsigned int    indx;
};

template<class outstream,class data_container_type>
outstream& operator<<(outstream& out,_ising_iterator<data_container_type>& it) {
  out<<it.index()<<" "<<*it;
  return out;
}
struct isingLattice1D {

  typedef int                                      value_type;
  typedef std::vector<value_type>                  data_container_type;
  typedef _ising_iterator<data_container_type>     iterator;

  template<class random_type>
  isingLattice1D(int L,random_type& ran) : data(L,1) {
	iterator it;
	for(it=begin();it!=end();++it){
	  double r = (ran()*1.0 / ran.max());
	  if(r > 0.5)
		(*it) = 1;
	  else
		(*it) = -1;
	  }
	calculateMagnetization();
	}
  template<class random_type>
  iterator random_site(random_type& ran) {
    
    unsigned int loc = static_cast<unsigned int>
                         (data.size()*(ran()*1.0/ran.max()));
    return iterator(data.begin()+loc,loc);
  }
  int effective_field(iterator iit) {
    int ret=0;
	iterator last = end(); --last;
	iterator current = begin(); current = iit;
	if(iit==begin()){
		iterator next = iit;
		++next;
		ret += *next;
		ret += *last;
	}
	else if (iit==last){
		iterator previous = iit;
		--previous;
		ret += *previous;
		ret += *(begin());
	}
	else {
		iterator next = iit;
		++next;
		ret += *next;
		iterator previous = iit;
		--previous;
		ret += *previous;
	}
	return ret;
  }
  
  void flip(iterator it){
	(*it) = -(*it);
	mag += 2*(*it);
	}

  iterator begin() { return iterator(data.begin(),0); }
  iterator end()   { return iterator(data.end(),data.size()); }
  unsigned int size() const { return data.size(); }
    
  double calculateMagnetization(){
	double s = std::accumulate(begin(),end(),0.);
	this->mag = s;
    return s;
	}
  
  double magnetization() { 
    return mag;
  }
  
  double magSign(){
	return (mag < 0) ? -1 : (mag > 0);
	}

  private:
    data_container_type	data;
	double mag;
};


#endif

