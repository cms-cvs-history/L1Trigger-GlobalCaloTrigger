#ifndef GCTTEST_H_
#define GCTTEST_H_

/*!
 * \class gctTestFunctions
 * \brief Top level skeleton for standalone testing
 * 
 * Contains a bunch of test functionality migrated from 
 * standalone test programs, to reduce duplication of code
 * and allow the tests to be done in an EDAnalyzer.
 *
 * Split into separate files to keep file sizes down.
 *
 * \author Greg Heath
 * \date March 2007
 *
 */

#include <vector>

class L1CaloEmCand;
class L1CaloRegion;

class gctTestElectrons;
class gctTestEnergyAlgos;
class gctTestFirmware;
class gctTestHtAndJetCounts;
class gctTestHfEtSums;

class L1GlobalCaloTrigger;

class gctTestFunctions
{
public:

  // structs and typedefs

  // Constructor and destructor
  gctTestFunctions();
  ~gctTestFunctions();

  /// Clear vectors of input data
  void reset();

  /// Load another event into the gct. Overloaded for the various ways of doing this.
  void loadNextEvent(L1GlobalCaloTrigger* &gct, const bool simpleEvent, const int16_t bx);
  void loadNextEvent(L1GlobalCaloTrigger* &gct, const std::string fileName, bool &endOfFile, const int16_t bx);
  void loadNextEvent(L1GlobalCaloTrigger* &gct, const std::string fileName, const int16_t bx);

  /// Read the input electron data (after GCT processing).
  void fillElectronData(const L1GlobalCaloTrigger* gct);

  /// Read the firmware results from a file for the next event
  void fillJetsFromFirmware(const std::string &fileName);

  /// Read the input jet data from the jetfinders (after GCT processing).
  void fillRawJetData(const L1GlobalCaloTrigger* gct);

  /// Check the electron sorter
  bool checkElectrons(const L1GlobalCaloTrigger* gct) const;

  /// Check the jet finder against results from the firmware
  bool checkJetFinder(const L1GlobalCaloTrigger* gct) const;

  /// Check the energy sums algorithms
  bool checkEnergySums(const L1GlobalCaloTrigger* gct) const;

  /// Check the Ht summing algorithms
  bool checkHtSums(const L1GlobalCaloTrigger* gct) const;

  /// Check the jet counting algorithms
  bool checkJetCounts(const L1GlobalCaloTrigger* gct) const;

  /// Check the Hf Et sums
  bool checkHfEtSums(const L1GlobalCaloTrigger* gct) const;

  /// Analyse calculation of energy sums in firmware
  bool checkEnergySumsFromFirmware(const L1GlobalCaloTrigger* gct, const std::string &fileName) const;

private:

  gctTestElectrons*      theElectronsTester;
  gctTestEnergyAlgos*    theEnergyAlgosTester;
  gctTestFirmware*       theFirmwareTester;
  gctTestHtAndJetCounts* theHtAndJetCountsTester;
  gctTestHfEtSums*       theHfEtSumsTester;

  std::vector< std::vector<L1CaloEmCand> > m_inputEmCands;
  std::vector< std::vector<L1CaloRegion> > m_inputRegions;

  int m_bxStart;
  int m_numOfBx;

  void bxRangeUpdate(const int16_t bx);

};

#endif /*GCTTEST_H_*/
