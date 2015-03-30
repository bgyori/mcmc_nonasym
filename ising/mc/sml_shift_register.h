#ifndef SHIFT_REGISTER_H
#define SHIFT_REGISTER_H

/*
 * sml headers
 */
#include <sml_bitwise_operations.h>
#include <sml_smart_array.h>

/**
  * algorithm to generate random sequence using
  *
  * random[i] = random[i-p] ^ random[i-q]
  *
  * for some values of p and q and with some initializations */


#include <limits>
#include <iostream>
#include <stdexcept>


namespace sml { /* namespace for standard mathematica library */


// forward declarations 

template<unsigned int,unsigned int,unsigned int> struct shift_register_default_seed;
template<unsigned int,unsigned int>              struct error_pq_values;
template<unsigned int>                           struct val_type;


/** 
  * template parameters 
  */
// seeder is a template for making the initializations, note the
// use of template template parameter to enforce consistent types.
// p and q are variables for shift register, possible values for
// p and q are:
// ---------------------------------------------
//    p          q
// ---------------------------------------------
//   89         38
//  127         1, 7, 15 30, 63
//  521         32, 48, 158, 168
//  607         105, 147, 273
// 1279         216, 418
// 2281         715, 915, 1029
// 3217         67, 576
// 4423         271, 369, 370, 649, 1393, 1419
// 9689         84, 471, 1836, 2444, 4187
// ---------------------------------------------


template<
	 unsigned int p_value,
	 unsigned int q_value,
	 unsigned int nbit = 8*sizeof(unsigned long),

         template<unsigned int, unsigned int, unsigned int> 
	 class seeder = shift_register_default_seed
        >

struct shift_register { 

    // determine the value type according to nbit
    typedef typename val_type<nbit>::result  value_type;
    
    // make template parameters avialable to public
    enum { nbits = nbit, p = p_value, q = q_value };

    // -----------------------------------------------------
    // constructor
    shift_register(const value_type init) : ran_vec(vec_size), seed_vec(p), iy(q) { 

      seeder<p,q,nbit> seed;
      seed(init,seed_vec);
      counter = 0;
      generate();
      // max value given by "000111...11" nbits of `1's.
      max_value = 0;
      value_type one = 1;
      for(value_type i=0;i<nbit;++i) {
	max_value = bit_or(bit_left_shift(  max_value, one ), one);
      }
    }

    // -----------------------------------------------------
    const value_type& operator()() {
      if(counter==vec_size) { counter=0; generate(); }
      return ran_vec[(counter++)];
    }

    // -----------------------------------------------------
    value_type min() const { return 0; }
    value_type max() const { return max_value; }

    // -----------------------------------------------------
    void test();

  private:

    // size of array for ran_vec 
    enum { vec_size = q*256 };

    // to detect errors in p and q values, if p and q
    // values are wrong, this line will give compile time
    // error
    error_pq_values<p,q>  undefined_pq_values;

    //value_type    seed_vec[p],ran_vec[vec_size];
    typename smart_array<value_type>::type ran_vec;
    typename smart_array<value_type>::type seed_vec;

    // used in generate function
    typename smart_array<value_type>::type iy;

    unsigned int  counter;
    value_type    max_value;
   
    // -----------------------------------------------------
    // use to generate the random sequence of size vec_size

    void generate() {

        enum { r = p - q, repeat = vec_size/q };

	//int iy[q];

	for(value_type j=0; j<repeat; j++){
	  for(value_type i=0; i<q; i++){
	    iy[i] = bit_xor(seed_vec[i], seed_vec[i+r]);
	    ran_vec[j*q+i]=iy[i];
	  }

	  for(value_type i=0; i<r; i++){
	    seed_vec[i]=seed_vec[i+q];
	  }
	  for(value_type i=0; i<q; i++){
	    seed_vec[i+r]=iy[i];
	  }
	}
    }
    // ---------------------------------------------------
};  /* end of class shift_register */



// ----------------------------------------------------
/**
  * the default algorithm to generate seed for shift register
  * pesudo-random number generator
  */
// ----------------------------------------------------


template<
	 unsigned int p,
	 unsigned int q,
	 unsigned int nbit
        >
struct shift_register_default_seed { 

    typedef typename val_type<nbit>::result  value_type;

    shift_register_default_seed ( ) { 
      mask = 0;
      value_type one = 1;
      for(value_type i=0;i<nbit;++i) {
	mask = bit_or(bit_left_shift( mask, one ), one);
      }
    }
    // ----------------------------------------------------
    // irseed must have memory size of "p" allocated

    template<class container>
    void operator()(const value_type& init,container& irseed) const {

        const value_type one = 1;
	value_type pp,qq;
	typename smart_array<value_type>::type iy(p);

	pp = p - 2;
	qq = p - 2+ q;

	for(value_type i=0;  i<p;   i++) { irseed[i] = 0;                                 }
	for(value_type i=0;  i<=30; i++) { iy[i] = bit_and(bit_right_shift(init, i), one);}
	for(value_type i=31; i<p;   i++) { iy[i] = bit_xor(iy[i-31], iy[i-13]);           }

	for(value_type i=0; i<=p*nbit-1; i++){
	  pp = (pp+1) % (p-1);
	  qq = (qq+1) % (p-1);
	  iy[pp] = bit_xor(iy[pp], iy[qq]);
	}
	for(value_type i=0; i<p; i++){
	  for(value_type j=0; j<=nbit-1; j++){
	    pp = (pp+1) % (p-1);
	    qq = (qq+1) % (p-1);
	    iy[pp] = bit_xor(iy[pp], iy[qq]);
	    irseed[i] = bit_or(bit_left_shift(irseed[i], one), iy[pp]);
	  }
	}
	pp = p - 2;
	qq = p - 2+ q;

	for(value_type i=1; i<=100000; i++){
	    pp = (pp+1) % (p-1);
	    qq = (qq+1) % (p-1);
	    irseed[pp] = bit_xor(irseed[pp], irseed[qq]);
	}

        for(value_type i=0;i<p;i++) { irseed[i] = bit_and(irseed[i], mask); }
    }
    // --------------------------------

    void test() {

      // mainly to test for right value of mask
      value_type sum=0, the_mask = mask;

      value_type one = 1;
      for(value_type i=0;i<8*sizeof(value_type);++i) {
        sum     += bit_and(the_mask,  one);
	the_mask = bit_right_shift(the_mask, one);
      }

      try {
	if(sum!=nbit) {
	  throw std::logic_error("seed test failed on mask");
	}
      }
      catch (std::logic_error le) {
	std::cerr<<__FILE__<<__LINE__<<" "<<le.what()<<std::endl;
	throw;
      }

    }
    // --------------------------------

  private:

    value_type mask;

};  /* end of class shift_register_default_seed */


/* ************************************

   determine the value_type base on nbit

   ************************************ */

// --------------------------------------------------------
// a class for checking if 64 bit is avialiable in this machine
template<bool> struct error_check64;
template<>     struct error_check64<true> {};

// --------------------------------------------------------
// returns unsigned long if template arguement is true
// otherwise do a compile time check and returns unsigned long long

template<bool> struct val32or64;

// value_type for 32 bits and less
template<>
struct val32or64<true>  { typedef unsigned int  result; };

// value_type for 64 bits and less
template<>
struct val32or64<false> { 

  typedef unsigned long long val_t;

  enum { max_size = sizeof(val_t)*8 };

  typedef error_check64<false> no_false_bit_available;

  typedef error_check64<max_size==64> no_64_bit_available;
  typedef val_t                       result; 
};
// --------------------------------------------------------

// val_type<nbit>::results give unsigned int if nbit is less than 32
// and give unsigned long long if nbit > 32 and less than or equal 64.
// gives compile time error if nbit > 64.

template<unsigned int nbit>
struct val_type {

  // make sure nbit is less or equal to 64
  error_check64<(nbit<=64)>    nbit_more_than_64_error;

  // define the required value type
  typedef typename val32or64<(nbit<=32)>::result  result;
};
// --------------------------------------------------------

/* ************************************

   this code enables compile time check
   for errors in p and q values

   ************************************ */
template<unsigned int p,unsigned int q> struct error_pq_values;
template<> struct error_pq_values<  89,  38> { };
template<> struct error_pq_values< 127,   1> { };
template<> struct error_pq_values< 127,   7> { };
template<> struct error_pq_values< 127,  15> { };
template<> struct error_pq_values< 127,  30> { };
template<> struct error_pq_values< 127,  63> { };
template<> struct error_pq_values< 521,  32> { };
template<> struct error_pq_values< 521,  48> { };
template<> struct error_pq_values< 521, 158> { };
template<> struct error_pq_values< 521, 168> { };
template<> struct error_pq_values< 607, 147> { };
template<> struct error_pq_values< 607, 273> { };
template<> struct error_pq_values<1279, 216> { };
template<> struct error_pq_values<1279, 418> { };
template<> struct error_pq_values<2281, 715> { };
template<> struct error_pq_values<2281, 915> { };
template<> struct error_pq_values<2281,1029> { };
template<> struct error_pq_values<3217,  67> { };
template<> struct error_pq_values<3217, 576> { };
template<> struct error_pq_values<4423, 271> { };
template<> struct error_pq_values<4423, 369> { };
template<> struct error_pq_values<4423, 370> { };
template<> struct error_pq_values<4423, 649> { };
template<> struct error_pq_values<4423,1393> { };
template<> struct error_pq_values<4423,1419> { };
template<> struct error_pq_values<9689,  84> { };
template<> struct error_pq_values<9689, 471> { };
template<> struct error_pq_values<9689,1836> { };
template<> struct error_pq_values<9689,2444> { };
template<> struct error_pq_values<9689,4187> { };

typedef shift_register<9689,4187>  shift_register_default;

// ----------------------------------------------------------
// the test code for shift_register.h

template<
	 unsigned int p,
	 unsigned int q,
	 unsigned int nbit,
         template<unsigned int, unsigned int, unsigned int> 
	 class seeder
        >
void
shift_register<p,q,nbit,seeder>::test() {

      // test the seed
      seeder<p,q,nbit>  seed;
      try {
	seed.test();
      }
      catch (...) {
	std::cerr<<"test failed from "<<__FILE__<<__LINE__<<std::endl;
	throw;
      }

      // first store the current state of the system so that
      // test will not affect results
      const unsigned int stored_counter = counter;
      typename smart_array<value_type>::type stored_ran_vec(ran_vec);
      typename smart_array<value_type>::type stored_seed_vec(seed_vec);


      // generate a new sequence
      generate();

      // flag for testing
      bool failed = false;

      // check bounds
      for(value_type i=0;i<vec_size;++i) {
	try {
	  if(ran_vec[i]>max() || ran_vec[i]<min()) {
	    throw std::range_error("test failed out-of-bounds");
	  }
	} catch(std::range_error re) {
	  std::cerr<<re.what()<<" on i="<<i
	           <<" ran_vec[i]="<<ran_vec[i]
		   <<" range=["<<min()<<","<<max()<<"]"<<std::endl;
	  failed = true;
	}
      }

      // algorithm is random[i] = random[i-p] ^ random[i-q]
      // choose a few values and test this xor
      // operation, choose 2*p, 3*p/2, 0, vec_size-1
      enum { r = p-q };

      try {
        if(ran_vec[2*p]!=bit_xor(ran_vec[p],ran_vec[2*p-q])) {
	  throw std::logic_error("test failed on i=2p");
        }
      } catch(std::logic_error le) {
	std::cerr<<le.what()<<std::endl;
	failed = true;
      }

      try {
	if(ran_vec[3*p/2]!=bit_xor(ran_vec[(3*p/2)-p],ran_vec[(3*p/2)-q])) {
	  throw std::logic_error("test failed on i=3p/2");
	}
      } catch(std::logic_error le) {
	std::cerr<<le.what()<<std::endl;
	failed = true;
      }

      try {
	if(ran_vec[0]!=bit_xor(stored_ran_vec[vec_size-p],stored_ran_vec[vec_size-q])) {
	  throw std::logic_error("test failed on i=0");
	}
      } catch(std::logic_error le) {
	std::cerr<<le.what()<<std::endl;
	failed = true;
      }

      try {
	if(ran_vec[vec_size-1]!=bit_xor(ran_vec[vec_size-1-p],ran_vec[vec_size-1-q])) {
	  throw std::logic_error("test failed on i=vec_size-1");
	}
      } catch(std::logic_error le) {
	std::cerr<<le.what()<<std::endl;
	failed = true;
      }

      // restore to state after test
      counter = stored_counter;
      for(value_type i=0;i<stored_ran_vec.size();++i) {
	ran_vec[i] = stored_ran_vec[i];
      }
      for(value_type i=0;i<stored_seed_vec.size();++i) {
	seed_vec[i] = stored_seed_vec[i];
      }

      if(failed) throw std::logic_error("shift_register test failed");
      else std::cerr<<__FILE__<<__LINE__
                    <<" shift_register.h test passed "<<std::endl;

    }



} /* namespace for standard mathematica library */ 

#endif

