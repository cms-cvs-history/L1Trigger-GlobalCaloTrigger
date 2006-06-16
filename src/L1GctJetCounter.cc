#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCounter.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinder.h"
 
#include "FWCore/Utilities/interface/Exception.h"  

#include <iostream>
using namespace std;

//DEFINE STATICS
const unsigned int L1GctJetCounter::MAX_LEAF_CARDS = L1GctWheelJetFpga::MAX_LEAF_CARDS;
const unsigned int L1GctJetCounter::MAX_JETS_PER_LEAF = L1GctJetLeafCard::MAX_JET_FINDERS*
                                                        L1GctJetFinder::MAX_JETS_OUT;
const unsigned int L1GctJetCounter::MAX_JETS_TO_COUNT = L1GctJetCounter::MAX_LEAF_CARDS*
                                                        L1GctJetCounter::MAX_JETS_PER_LEAF;


L1GctJetCounter::L1GctJetCounter(int id, vector<L1GctJetLeafCard*> leafCards,
                               L1GctJetCounterLut* jetCounterLut):
  m_id(id),
  m_jetLeafCards(leafCards),
  m_jetCounterLut(jetCounterLut),
  m_jets(MAX_JETS_TO_COUNT)
{
  //Check jetfinder setup
  if(m_id < 0 || m_id%100 >= 12 || m_id/100 >= 2)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetCounter::L1GctJetCounter() : Jet Counter ID " << m_id << " has been incorrectly constructed!\n"
    << "ID number should be between the range of 0 to 11, or 100 to 111\n";
  } 
  
  if(m_jetLeafCards.size() != MAX_LEAF_CARDS)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetCounter::L1GctJetCounter() : Jet Counter ID " << m_id << " has been incorrectly constructed!\n"
    << "This class needs " << MAX_LEAF_CARDS << " leaf card pointers, yet only " << m_jetLeafCards.size()
    << " leaf card pointers are present.\n";
  }
  
  for(unsigned int i = 0; i < m_jetLeafCards.size(); ++i)
  {
    if(m_jetLeafCards[i] == 0)
    {
      throw cms::Exception("L1GctSetupError")
      << "L1GctJetCounter::L1GctJetCounter() : Jet Counter ID " << m_id << " has been incorrectly constructed!\n"
      << "Leaf card pointer " << i << " has not been set!\n";
    }
  }
  
  if(m_jetCounterLut == 0)
  {
    throw cms::Exception("L1GctSetupError")
    << "L1GctJetCounter::L1GctJetCounter() : Jet Counter ID " << m_id << " has been incorrectly constructed!\n"
    << "The jet counter LUT pointer has not been set!\n";  
  }
}

L1GctJetCounter::~L1GctJetCounter()
{
}

ostream& operator << (ostream& os, const L1GctJetCounter& algo)
{
  os << "===L1GctJetCounter===" << endl;
  os << "ID = " << algo.m_id << endl;
  os << "JetCounterLut* = " <<  algo.m_jetCounterLut << endl;
  os << *algo.m_jetCounterLut << endl;
  os << "No of Leaf cards " << algo.m_jetLeafCards.size() << endl;
  for (unsigned i=0; i<algo.m_jetLeafCards.size(); i++) {
    os << "JetLeafCard* " << i << " = " << algo.m_jetLeafCards[i]<< endl;
    os << "No of jets from this leafCard, jetFinder A " << algo.m_jetLeafCards.at(i)->getOutputJetsA().size() << endl;
    os << "No of jets from this leafCard, jetFinder B " << algo.m_jetLeafCards.at(i)->getOutputJetsB().size() << endl;
    os << "No of jets from this leafCard, jetFinder C " << algo.m_jetLeafCards.at(i)->getOutputJetsC().size() << endl;
  }
  os << "No of input jets " << algo.m_jets.size() << endl;
//   for(unsigned i=0; i < algo.m_jets.size(); ++i)
//     {
//       os << algo.m_jets[i]; 
//     }
  os << "Value of counter " << algo.m_value << endl;
  os << endl;

  return os;
}


void L1GctJetCounter::reset()
{
  for (unsigned i=0; i<m_jets.size(); i++) {
    m_jets.at(i).setupJet(0, 0, 0, true);
  }
  m_value.reset();
}

// Load the m_jets vector
void L1GctJetCounter::fetchInput()
{
  int jetnum=0;
  for (unsigned i=0; i<m_jetLeafCards.size(); i++) {
    if (jetnum+MAX_JETS_PER_LEAF>=m_jets.size()) {
      throw cms::Exception("L1GctProcessingError")
	<< "L1GctJetCounter id= " << m_id << " trying to input too many jets for Leaf Card number " << i << endl;
    }
    L1GctJetLeafCard* jlc = m_jetLeafCards.at(i);
    for (int j=0; j<L1GctJetFinder::MAX_JETS_OUT; j++) {
      m_jets.at(jetnum++) = jlc->getOutputJetsA().at(j);
      m_jets.at(jetnum++) = jlc->getOutputJetsB().at(j);
      m_jets.at(jetnum++) = jlc->getOutputJetsC().at(j);
    }
  }
}

/// set the m_jets vector for test purposes
void L1GctJetCounter::setJets(JetVector jets)
{
  m_jets = jets;
}

/// Count the jets passing cuts
void L1GctJetCounter::process() 
{
  for (unsigned i=0; i<m_jets.size(); i++) {
    if (m_jetCounterLut->passesCut(m_jets.at(i))) { m_value++; }
  }
}