#include "L1Trigger/GlobalCaloTrigger/interface/L1GctHfEtSumsLut.h"

//DEFINE STATICS
const int L1GctHfEtSumsLut::NAddress=L1GctHfLutSetup::kHfEtSumBits;
const int L1GctHfEtSumsLut::NData   =L1GctHfLutSetup::kHfOutputBits;

L1GctHfEtSumsLut::L1GctHfEtSumsLut(const L1GctHfLutSetup::hfLutType& type, const L1GctHfLutSetup* const fn) :
  L1GctLut<NAddress,NData>(),
  m_lutFunction(fn),
  m_lutType(type)
{
  if (fn != 0) m_setupOk = true;
}

L1GctHfEtSumsLut::L1GctHfEtSumsLut(const L1GctHfLutSetup::hfLutType& type) :
  L1GctLut<NAddress,NData>(),
  m_lutFunction(0),
  m_lutType(type)
{
}

L1GctHfEtSumsLut::L1GctHfEtSumsLut() :
  L1GctLut<NAddress,NData>(),
  m_lutFunction(0),
  m_lutType()
{
}

L1GctHfEtSumsLut::L1GctHfEtSumsLut(const L1GctHfEtSumsLut& lut) :
  L1GctLut<NAddress,NData>(),
  m_lutFunction(lut.lutFunction()),
  m_lutType(lut.lutType())
{
}

L1GctHfEtSumsLut::~L1GctHfEtSumsLut()
{
}

L1GctHfEtSumsLut L1GctHfEtSumsLut::operator= (const L1GctHfEtSumsLut& lut)
{
  L1GctHfEtSumsLut temp(lut);
  return temp;
}

std::ostream& operator << (std::ostream& os, const L1GctHfEtSumsLut& lut)
{
  os << "===L1GctHfEtSumsLut===" << std::endl;
  os << "\n===Lookup table contents===\n" << std::endl;
  const L1GctLut<L1GctHfEtSumsLut::NAddress,L1GctHfEtSumsLut::NData>* temp=&lut;
  os << *temp;
  return os;
}

template class L1GctLut<L1GctHfEtSumsLut::NAddress,L1GctHfEtSumsLut::NData>;


uint16_t L1GctHfEtSumsLut::value (const uint16_t lutAddress) const
{
  return m_lutFunction->outputValue(m_lutType, lutAddress) ;
}
