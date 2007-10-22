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
  typedef std::vector<L1CaloRegion> RegionsVector;
  typedef std::vector<L1GctJetCand> JetsVector;
  typedef std::vector<L1GctJet>     RawJetsVector;

  struct rawJetData { RawJetsVector jets; unsigned htSum; };

  // Constructor and destructor
  gctTestHtAndJetCounts();
  ~gctTestHtAndJetCounts();

  /// Read the input jet data from the jetfinders (after GCT processing).
  void fillRawJetData(const L1GlobalCaloTrigger* gct);

  /// Check the Ht summing algorithms
  bool checkHtSums(const L1GlobalCaloTrigger* gct) const;

  /// Check the jet counting algorithms
  bool checkJetCounts(const L1GlobalCaloTrigger* gct) const;

private:

  //
  // FUNCTION PROTOTYPES FOR HT SUM CHECKING
  /// Works out the Ht for a jet from the raw input regions
  unsigned jetHtSum(const L1GctJetFinderBase* jf, const int jn) const;
  //=========================================================================

  //
  // FUNCTION PROTOTYPE FOR JET COUNTING
  /// Counts jets in cuts
  unsigned countJetsInCut(const std::vector<rawJetData>& jetList,
                          const L1GctJetCounterSetup::cutsListForJetCounter& cutList,
                          const L1GctJetEtCalibrationLut* lut) const;
  //=========================================================================

  rawJetData rawJetFinderOutput(const L1GctJetFinderBase* jf) const;

  std::vector<rawJetData> minusWheelJetDta;
  std::vector<rawJetData> plusWheelJetData;

};

#endif /*GCTTEST_H_*/
