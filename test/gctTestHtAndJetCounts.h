#ifndef GCTTESTHTANDJETCOUNTS_H_
#define GCTTESTHTANDJETCOUNTS_H_

/*!
 * \class gctTestHtAndJetCounts
 * \brief Test of the Ht and jet counts
 * 
 * Ht and jet count sum test functionality migrated from standalone test programs
 *
 * \author Greg Heath
 * \date March 2007
 *
 */
 
#include "CondFormats/L1TObjects/interface/L1GctJetCounterSetup.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctJetCand.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJet.h"

#include <vector>

class L1GlobalCaloTrigger;
class L1GctJetLeafCard;
class L1GctJetFinderBase;

class gctTestHtAndJetCounts
{
public:

  // structs and typedefs
  typedef std::vector<L1GctJetCand> JetsVector;
  typedef std::vector<L1GctJet>     RawJetsVector;

  struct rawJetData {
    RawJetsVector jets; unsigned htStripSum0; unsigned htStripSum1; bool htOverFlow;

    rawJetData() :
      jets(), htStripSum0(0), htStripSum1(0), htOverFlow(false) {}
    rawJetData(const RawJetsVector jv, const unsigned ht0, const unsigned ht1, const bool of) : 
      jets(jv), htStripSum0(ht0), htStripSum1(ht1), htOverFlow(of) {}
  };

  // Constructor and destructor
  gctTestHtAndJetCounts();
  ~gctTestHtAndJetCounts();

  /// Set array sizes for the number of bunch crossings
  void setBxRange(const int bxStart, const int numOfBx);

  /// Read the input jet data from the jetfinders (after GCT processing).
  void fillRawJetData(const L1GlobalCaloTrigger* gct);

  /// Check the Ht summing algorithms
  bool checkHtSums(const L1GlobalCaloTrigger* gct) const;

  /// Check the jet counting algorithms
  bool checkJetCounts(const L1GlobalCaloTrigger* gct) const;

private:

  //

  //
  // FUNCTION PROTOTYPE FOR JET COUNTING
  /// Counts jets in cuts
  unsigned countJetsInCut(const std::vector<rawJetData>& jetList,
                          const L1GctJetCounterSetup::cutsListForJetCounter& cutList,
                          const L1GctJetEtCalibrationLut* lut,
			  const int bx) const;
  //=========================================================================

  rawJetData rawJetFinderOutput(const L1GctJetFinderBase* jf, const int bx) const;

  int m_bxStart;
  int m_numOfBx;

  std::vector<rawJetData> minusWheelJetDta;
  std::vector<rawJetData> plusWheelJetData;

  int etComponent(const unsigned Emag0, const unsigned fact0,
		  const unsigned Emag1, const unsigned fact1) const;
};

#endif /*GCTTEST_H_*/
