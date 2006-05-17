#ifndef L1GCTETTYPES_H
#define L1GCTETTYPES_H

#include <boost/cstdint.hpp>
#include <ostream>


/* Signed integer representations */

// class to store a 2s complement number
// in the range -X to X-1
// including overflow

template <int nBits>
class L1GctTwosComplement {
 public:
  L1GctTwosComplement();
  L1GctTwosComplement(uint32_t raw);
  L1GctTwosComplement(int value);
  ~L1GctTwosComplement() { }

  // reset value and overflow to zero
  void reset();

  // set the raw data
  void setRaw(uint32_t raw);

  // set value from signed int
  void setValue(int value);

  // set the overflow bit
  void setOverFlow(bool oflow) { m_overFlow = oflow; }

  // access raw data
  uint32_t raw() const { return m_data; }

  // access value as signed int
  int value() const;

  // access overflow
  bool overFlow() const { return m_overFlow; }

  // return number of bits
  int size() const { return m_nBits; }

  // add two numbers of the same size
  L1GctTwosComplement operator+ (const L1GctTwosComplement &rhs) const;

  // overload = operator
  L1GctTwosComplement& operator= (int value);

 private:

  // number of bits (for overflow checking)
  int m_nBits;

  // the raw data
  uint32_t m_data;

  // the overflow bit
  bool m_overFlow;

  static const int MAX_NBITS = 24;
  static const int MAX_VALUE = 1<<(MAX_NBITS-1);

  // function to check overflow
  void checkOverFlow(uint32_t rawValue, uint32_t &maskValue, bool &overFlow);
};

// construct with # bits and set to zero
template <int nBits>
L1GctTwosComplement<nBits>::L1GctTwosComplement() {
  m_nBits = nBits>0 && nBits<MAX_NBITS ? nBits : 16 ;
  this->reset();
}

// construct from # bits and raw data 
template <int nBits>
L1GctTwosComplement<nBits>::L1GctTwosComplement(uint32_t raw) {
  m_nBits = nBits>0 && nBits<MAX_NBITS ? nBits : 16 ;
  this->setRaw(raw);
}

// construct from # bits and value
template <int nBits>
L1GctTwosComplement<nBits>::L1GctTwosComplement(int value) {
  m_nBits = nBits>0 && nBits<MAX_NBITS ? nBits : 16 ;
  this->setValue(value);
}

// construct from reset value and overflow to zero
template <int nBits>
void L1GctTwosComplement<nBits>::reset() {
  m_data = static_cast<uint32_t>(0);
  m_overFlow = false;
}

// set value from uint32_t
template <int nBits>
void L1GctTwosComplement<nBits>::setRaw(uint32_t raw) {
  checkOverFlow(raw, m_data, m_overFlow);
}

// set value from int
template <int nBits>
void L1GctTwosComplement<nBits>::setValue(int value) {
  int chkValue, posValue;
  uint32_t raw;

  // Make sure we have an integer in the range MAX_NBITS
  chkValue = value;
  if (chkValue<-MAX_VALUE) { chkValue =  -MAX_VALUE; m_overFlow = true; }
  if (chkValue>=MAX_VALUE) { chkValue = MAX_VALUE-1; m_overFlow = true; }

  // Transform negative values to large positive values
  posValue = chkValue<0 ? chkValue + (1<<MAX_NBITS) : chkValue ;
  raw = static_cast<uint32_t>(posValue);

  // Use the setRaw method to check overflow for our given size nBits
  this->setRaw(raw);
}

// return value as int
template <int nBits>
int L1GctTwosComplement<nBits>::value() const {
  int value, result;
  int maxValueInNbits;
  maxValueInNbits = 1<<(m_nBits-1);
  value  = static_cast<int>(m_data);
  result = value < maxValueInNbits ? value : value - (1<<MAX_NBITS) ;
  return result;
}

// add two numbers
template <int nBits>
L1GctTwosComplement<nBits>
L1GctTwosComplement<nBits>::operator+ (const L1GctTwosComplement<nBits> &rhs) const {

  // temporary variable for storing the result (need to set its size)
  L1GctTwosComplement<nBits> temp;
  uint32_t sum;
  bool ofl;

  // do the addition here
  sum = this->raw() + rhs.raw();
  ofl = this->overFlow() || rhs.overFlow();

  //fill the temporary argument
  temp.setRaw(sum);
  temp.setOverFlow(temp.overFlow() || ofl);

  // return the temporary
  return temp;

}

// overload assignment by int
template <int nBits>
L1GctTwosComplement<nBits>& L1GctTwosComplement<nBits>::operator= (int value) {
  
  this->setValue(value);
  return *this;

}

// Here's the check overflow function
template <int nBits> 
void L1GctTwosComplement<nBits>::checkOverFlow(uint32_t rawValue, uint32_t &maskValue, bool &overFlow) {
  uint32_t signBit = 1<<(m_nBits-1);
  uint32_t signExtendBits = (static_cast<uint32_t>(MAX_VALUE)-signBit)<<1;
  uint32_t value;
  bool ofl;

  if ((rawValue&signBit)==0) {
    value = rawValue & ~signExtendBits;
  } else {
    value = rawValue | signExtendBits;
  }
  ofl = value != rawValue;

  maskValue = value;
  overFlow  = ofl;

}


/* unsigned integer representations */

template <int nBits>
class L1GctUnsignedInt {

 public:

  L1GctUnsignedInt();
  L1GctUnsignedInt(unsigned value);
  ~L1GctUnsignedInt();

  // reset value and overflow to zero
  void reset() { m_value = static_cast<unsigned>(0); m_overFlow = false; }

  // Set value from unsigned
  void setValue(unsigned value);

  // set the overflow bit
  void setOverFlow(bool oflow) { m_overFlow = oflow; }

  // access value as unsigned
  unsigned value() const { return m_value; }

  // access overflow
  bool overFlow() const { return m_overFlow; }

  // return number of bits
  int size() const { return m_nBits; }

  // add two numbers
  L1GctUnsignedInt operator+ (const L1GctUnsignedInt &rhs) const;

 private:

  // number of bits
  int m_nBits;

  // value
  unsigned m_value;

  // overflow
  bool m_overFlow;

  static const int MAX_NBITS = 24;

};

/* unsigned integer representations */

template <int nBits>
L1GctUnsignedInt<nBits>::L1GctUnsignedInt() {

  m_nBits = nBits>0 && nBits<MAX_NBITS ? nBits : 16 ;
  this->reset();
}

template <int nBits>
L1GctUnsignedInt<nBits>::L1GctUnsignedInt(unsigned value) {

  m_nBits = nBits>0 && nBits<MAX_NBITS ? nBits : 16 ;
  this->setValue(value);
}

template <int nBits>
L1GctUnsignedInt<nBits>::~L1GctUnsignedInt()
{

}

// set value, checking for overflow
template <int nBits>
void L1GctUnsignedInt<nBits>::setValue(unsigned value)
{
  // check for overflow
  if (value >= (unsigned) 1<<(m_nBits-1) ) {
    m_overFlow = true;
  }

  // set value with bitmask
  m_value = value & (1<<(m_nBits-1) - 1);

}

// add two unsigneds
template <int nBits>
L1GctUnsignedInt<nBits>
L1GctUnsignedInt<nBits>::operator+ (const L1GctUnsignedInt<nBits> &rhs) const {

  // temporary variable for storing the result (need to set its size)
  L1GctUnsignedInt temp(nBits);

  unsigned sum;
  bool ofl;

  // do the addition here
  sum = this->value() + rhs.value();
  ofl = this->overFlow() || rhs.overFlow();

  //fill the temporary argument
  temp.setValue(sum);
  temp.setOverFlow(temp.overFlow() || ofl);

  // return the temporary
  return temp;

}

/// Here are the typedef's for the data types used in the energy sum calculations

typedef L1GctTwosComplement<12> L1GctEtComponent;
typedef L1GctUnsignedInt<12>    L1GctScalarEtVal;
typedef L1GctUnsignedInt<7>     L1GctEtAngleBin;


template <int nBits>
std::ostream& operator<<(std::ostream& s, const L1GctTwosComplement<nBits>& data);

template <int nBits>
std::ostream& operator<<(std::ostream& s, const L1GctUnsignedInt<nBits>& data);

// overload ostream<<
template <int nBits>
std::ostream& operator<<(std::ostream& s, const L1GctTwosComplement<nBits>& data) {

  s << "L1GctTwosComplement raw : " << data.raw() << ", " << "value : " << data.value();
  if (data.overFlow()) { s << " Overflow set! "; }
  s << std::endl;
  return s;

}

template <int nBits>
std::ostream& operator<<(std::ostream& s, const L1GctUnsignedInt<nBits>& data) {

  s << "L1GctUnsignedInt raw : " << data.raw() << ", " << "value : " << data.value();
  if (data.overFlow()) { s << " Overflow set! "; }
  s << std::endl;
  return s;

}

#endif