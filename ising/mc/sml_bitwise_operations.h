#ifndef BITWISE_OPERATION_H
#define BITWISE_OPERATION_H

namespace bitwise_op_type_check {

  template<class T1,class T2> struct is_same_type;
  template<class T>           struct is_same_type<T,T> {
    void operator()() { }
  };
}

/* namespace for standard mathematical library */
namespace sml {

template<class T1,class T2>
T1 bit_and(const T1& v1,const T2& v2) {
  bitwise_op_type_check::is_same_type<T1,T2> ist; ist();
  return ( v1 & v2 );
}

template<class T1,class T2>
T1 bit_or(const T1& v1,const T2& v2) {
  bitwise_op_type_check::is_same_type<T1,T2> ist; ist();
  return ( v1 | v2 );
}

template<class T1,class T2>
T1 bit_xor(const T1& v1,const T2& v2) {
  bitwise_op_type_check::is_same_type<T1,T2> ist; ist();
  return ( v1 ^ v2 );
}

template<class T1,class T2>
T1 bit_right_shift(const T1& v1,const T2& v2) {
  bitwise_op_type_check::is_same_type<T1,T2> ist; ist();
  return ( v1 >> v2 );
}

template<class T1,class T2>
T1 bit_left_shift(const T1& v1,const T2& v2) {
  bitwise_op_type_check::is_same_type<T1,T2> ist; ist();
  return ( v1 << v2 );
}

} /* namespace sml */

#endif


