#ifndef SMART_ARRAY_H
#define SMART_ARRAY_H

#include <vector>
#include <stdexcept>
#include <iostream>


/* *************************************
 * smart_array class checks for index out
 * of bounds if the macro DEBUG_VERSION is
 * defined. Otherwise smart_array is the 
 * same as std::vector. When DEBUG_VERSION
 * is defined, the program runs about 30%
 * slower.
 * *************************************/


/* a namespace for standard mathematical library */
namespace sml {

// forward declaration
template<class value_type,class unused_allocator> struct debug_array;


// a work around for template typedef

template<class value_type>
struct smart_array {

#if defined(DEBUG_VERSION)
  typedef debug_array<value_type,std::allocator<value_type> >  type;
#else
  typedef std::vector<value_type>  type;
#endif

};


// ---------------------------------------
template<class val_type,
         class allocator = std::allocator<val_type> >
struct debug_array : public std::vector<val_type,allocator> { 


  typedef std::vector<val_type,allocator> base_type;

  typedef typename base_type::value_type           value_type;
  typedef typename base_type::const_iterator       const_iterator;
  typedef typename base_type::iterator             iterator;

  debug_array()  { }

  debug_array(unsigned int s) : 
    std::vector<val_type,allocator>(s) { }

  debug_array(unsigned int s,val_type v) : 
    std::vector<val_type,allocator>(s,v) { }


  // --------------------------------------------------
  const val_type& operator[](unsigned int index) const {

    try {
      if(index<0 || index >= base_type::size()) {
	throw std::range_error("index out-of-range");
      }
    } catch (std::range_error re) {
      std::cerr<<re.what()<<" from "<<__FILE__<<__LINE__<<std::endl;
      throw;
    }

    return base_type::operator[](index);

  }

  // --------------------------------------------------
  val_type& operator[](unsigned int index) {

    try {
      if(index<0 || index >= base_type::size()) {
	throw std::range_error("index out-of-range");
      }
    } catch (std::range_error re) {
      std::cerr<<re.what()<<" from "<<__FILE__<<__LINE__<<std::endl;
      throw;
    }

    return base_type::operator[](index);
  }
};

} /* namespace sml */

#endif

