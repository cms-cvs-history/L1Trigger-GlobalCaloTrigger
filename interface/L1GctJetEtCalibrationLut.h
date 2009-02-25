#ifndef L1GCTJETETCALIBRATIONLUT_H_
#define L1GCTJETETCALIBRATIONLUT_H_

#define JET_ET_CAL_LUT_ADD_BITS 11
#define JET_ET_CAL_LUT_DAT_BITS 16

#include "L1Trigger/GlobalCaloTrigger/src/L1GctLut.h"

class L1GctJetFinderParams;
class L1CaloEtScale;

/*!
 * \author Robert Frazier & Greg Heath
 * \date May 2006
 */

/*! \class L1GctJetEtCalibrationLut
 * \brief Jet Et calibration LUT
 * 
 * Input is 10 bit Et and 4 bit eta
 * Outputs are 6 bit rank (for jet sorting) and 10 bit Et (for Ht calculation)
 * 
 * Modified March 2007 to remove the actual calculation to a separate class
 * Modified October 2008 to have separate LUTs for each eta, as in the firmware
 *
 */

class L1GctJetEtCalibrationLut : public L1GctLut<JET_ET_CAL_LUT_ADD_BITS,JET_ET_CAL_LUT_DAT_BITS>
{
 public:
  static const int NAddress;
  static const int NData;
  static const unsigned JET_ENERGY_BITWIDTH;   ///< Input bitwidth of jet energy; must be 10 or more

  L1GctJetEtCalibrationLut();
  virtual ~L1GctJetEtCalibrationLut();

  // set components
  void setFunction(const L1GctJetFinderParams * const lutfn);
  void setOutputEtScale(const L1CaloEtScale * const scale);
  void setEtaBin(const unsigned eta);

  // get components
  const L1GctJetFinderParams* getFunction() const { return m_lutFunction; }
  const L1CaloEtScale* getOutputEtScale() const { return m_outputEtScale; }
  unsigned etaBin() const { return static_cast<unsigned>(m_etaBin); }

  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctJetEtCalibrationLut& lut);
  
 protected:
  

  virtual uint16_t value (const uint16_t lutAddress) const;

 private:

  const L1GctJetFinderParams* m_lutFunction;
  const L1CaloEtScale * m_outputEtScale;

  uint8_t m_etaBin;

};

std::ostream& operator << (std::ostream& os, const L1GctJetEtCalibrationLut& lut);

#endif /*L1GCTJETETCALIBRATIONLUT_H_*/
