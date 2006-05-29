#ifndef L1GCTWHEELJETFPGA_H_
#define L1GCTWHEELJETFPGA_H_

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctProcessor.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetLeafCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEtTypes.h"

#include <vector>

class L1GctWheelJetFpga : public L1GctProcessor
{
public:
  typedef std::vector<L1GctJetCand> JetVector;
  static const int MAX_JETS_OUT;   ///< Max number of jets of each type we output.

  L1GctWheelJetFpga(int id, std::vector<L1GctJetLeafCard*> inputLeafCards);
  ~L1GctWheelJetFpga();

  /// Overload << operator
  friend std::ostream& operator << (std::ostream& os, const L1GctWheelJetFpga& fpga);

  /// clear internal buffers
  virtual void reset();

  /// get input data from sources
  virtual void fetchInput();

  /// process the data, fill output buffers
  virtual void process();

  /// set input data      
  void setInputJet(int i, L1GctJetCand jet); 
  void setInputHt (int i, unsigned ht);
    
  /// get the input jets. Jets 0-5 from leaf card 0, jetfinderA.  Jets 6-11 from leaf card 0, jetfinder B... etc.
  JetVector getInputJets() const { return m_inputJets; }
    
  /// get the input Ht
  L1GctScalarEtVal inputHt(unsigned leafnum) const { return m_inputHt[leafnum]; }
    
  /// get the output jets
  JetVector getCentralJets() const { return m_centralJets; }
  JetVector getForwardJets() const { return m_forwardJets; }
  JetVector getTauJets() const { return m_tauJets; }
    
  /// get the output Ht and jet counts
  L1GctScalarEtVal getOutputHt() const { return m_outputHt; }
  L1GctJcWheelType getOutputJc(unsigned jcnum) const { return m_outputJc[jcnum]; }

private:

  static const unsigned int MAX_LEAF_CARDS;  ///< Max number of leaf card pointers
  static const int MAX_JETS_IN;    ///< Maximum number of jets we can have as input
  static const int MAX_RAW_CJETS;  ///< Max. possible central jets to sort over
  static const int MAX_RAW_FJETS;  ///< Max. possible forward jets to sort over
  static const int MAX_RAW_TJETS;  ///< Max. possible tau jets to sort over
    
  /// algo ID
  int m_id;

  /// the jet leaf cards
  std::vector<L1GctJetLeafCard*> m_inputLeafCards;
    
  /// input data. Jets 0-5 from leaf card 0, jetfinderA.  Jets 6-11 from leaf card 0, jetfinder B... etc.
  JetVector m_inputJets;
    
  // Holds the all the various inputted jets, re-addressed using proper GCT->GT jet addressing
  JetVector m_rawCentralJets; 
  JetVector m_rawForwardJets; 
  JetVector m_rawTauJets; 

  // input Ht sums from each leaf card
  std::vector<L1GctScalarEtVal> m_inputHt;

  // output data
  JetVector m_centralJets;
  JetVector m_forwardJets;
  JetVector m_tauJets;
    
  // data sent to GlobalEnergyAlgos
  L1GctScalarEtVal m_outputHt;
  std::vector<L1GctJcWheelType> m_outputJc;
      
  //PRIVATE METHODS
  /// Puts the output from a jetfinder into the correct index range of the m_inputJets array. 
  void storeJets(JetVector jets, unsigned short iLeaf, unsigned short offset);
  /// Classifies jets into central, forward or tau, and re-addresses them using global co-ords.
  void classifyJets();
  /// Sizes the m_rawTauJetsVec, and then sets all tauVeto bits to false.
  void setupRawTauJetsVec();
};

std::ostream& operator << (std::ostream& os, const L1GctWheelJetFpga& fpga);

#endif /*L1GCTWHEELJETFPGA_H_*/
