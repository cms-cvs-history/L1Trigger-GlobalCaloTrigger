#ifndef L1GCTJETCOUNTER_H_
#define L1GCTJETCOUNTER_H_

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctProcessor.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetLeafCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEtTypes.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCounterLut.h"

#include <boost/cstdint.hpp> //for uint16_t
#include <vector>

/*!
 * \class L1GctJetCounter
 * \brief Counts jets in one Wheel that pass criteria encoded in a JetCounterLut
 *
 * The actual contents of this class are fairly simple, since
 * all the real work is done elsewhere.
 *  
 * \author Greg Heath
 * \date June 2006
 */



class L1GctJetCounter : public L1GctProcessor
{
public:
  //Typedefs
  typedef std::vector<L1GctJetCand> JetVector;

  //Statics
  static const unsigned int MAX_LEAF_CARDS;
  static const unsigned int MAX_JETS_PER_LEAF;
  static const unsigned int MAX_JETS_TO_COUNT;
    
  /// id needs to encode Wheel and jet count numbers
  L1GctJetCounter(int id, std::vector<L1GctJetLeafCard*> leafCards,
                  L1GctJetCounterLut* jetCounterLut);
                 
  ~L1GctJetCounter();
   
  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctJetCounter& algo);

  /// clear internal buffers
  virtual void reset();

  /// get input data from sources
  virtual void fetchInput();

  /// process the data, fill output buffers
  virtual void process();

  /// set the input jets (for test purposes)
  void setJets(JetVector jets);

  /// get the JetCounterLut
  L1GctJetCounterLut* getJetCounterLut() const { return m_jetCounterLut; }

  /// get the jets
  JetVector getJets() const { return m_jets; }

  /// get the value of the counter, for input into the jet count sums
  L1GctJcWheelType getValue() const { return m_value;}

private:

  //Statics

  /// algo ID
  int m_id;

  /// Jet Leaf Card pointers
  std::vector<L1GctJetLeafCard*> m_jetLeafCards;
	
  /// Jet Et Converstion LUT pointer
  L1GctJetCounterLut* m_jetCounterLut;

  /// The jets to be counted
  JetVector m_jets;
    
  /// The value of the counter
  L1GctJcWheelType m_value;

  //PRIVATE METHODS
  
};

std::ostream& operator << (std::ostream& os, const L1GctJetCounter& algo);

#endif /*L1GCTJETCOUNTER_H_*/