#include "L1Trigger/GlobalCaloTrigger/test/gctTestHtAndJetCounts.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinderBase.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCounter.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCounterLut.h"

#include <math.h>
#include <iostream>
#include <cassert>

using namespace std;

//=================================================================================================================
//
/// Constructor and destructor

gctTestHtAndJetCounts::gctTestHtAndJetCounts() {}
gctTestHtAndJetCounts::~gctTestHtAndJetCounts() {}

//=================================================================================================================
//
/// Read the input jet data from the jetfinders (after GCT processing).
void gctTestHtAndJetCounts::fillRawJetData(const L1GlobalCaloTrigger* gct) {

  minusWheelJetDta.clear();
  plusWheelJetData.clear();
  minusWheelJetDta.resize(9*m_numOfBx);
  plusWheelJetData.resize(9*m_numOfBx);

  int bx=m_bxStart;
  int mPos=0;
  int pPos=0;
  for (int ibx=0; ibx<m_numOfBx; ibx++) { 
    int leaf=0;

    // Minus Wheel
    for ( ; leaf<3; leaf++) {
      minusWheelJetDta.at(mPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderA(), bx);
      minusWheelJetDta.at(mPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderB(), bx);
      minusWheelJetDta.at(mPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderC(), bx);
    }

    // Plus Wheel
    for ( ; leaf<6; leaf++) {
      plusWheelJetData.at(pPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderA(), bx);
      plusWheelJetData.at(pPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderB(), bx);
      plusWheelJetData.at(pPos++) = rawJetFinderOutput(gct->getJetLeafCards().at(leaf)->getJetFinderC(), bx);
    }

    bx++;
  }

}

//=================================================================================================================
//
/// Set array sizes for the number of bunch crossings
void gctTestHtAndJetCounts::setBxRange(const int bxStart, const int numOfBx){
  // Allow the start of the bunch crossing range to be altered
  // without losing previously stored jet data
  rawJetData temp;
  for (int bx=bxStart; bx<m_bxStart; bx++) {
    minusWheelJetDta.insert(minusWheelJetDta.begin(), 9, temp);
    plusWheelJetData.insert(plusWheelJetData.begin(), 9, temp);
  }

  m_bxStart = bxStart;

  // Resize the vectors without clearing previously stored values
  minusWheelJetDta.resize(9*numOfBx);
  plusWheelJetData.resize(9*numOfBx);
  m_numOfBx = numOfBx;
}

//=================================================================================================================
//
/// Check the Ht summing algorithms
bool gctTestHtAndJetCounts::checkHtSums(const L1GlobalCaloTrigger* gct) const
{
  bool testPass = true;
  L1GctGlobalEnergyAlgos* myGlobalEnergy = gct->getEnergyFinalStage();

  std::vector<rawJetData>::const_iterator mJet=minusWheelJetDta.begin();
  std::vector<rawJetData>::const_iterator pJet=plusWheelJetData.begin();

  for (int bx=0; bx<m_numOfBx; bx++) {
    unsigned htMinusVl = 0;
    unsigned htPlusVal = 0;

    //
    // Check the Ht calculation (starting from the found jets)
    //--------------------------------------------------------------------------------------
    //
    // Minus Wheel
    int leaf=0;
    for ( ; leaf<3; leaf++) {
      unsigned leafHtSum = 0;
      for (int jf=0; jf<3; jf++) {
	assert (mJet != minusWheelJetDta.end());
	leafHtSum += mJet->htSum;
	mJet++;
      }
      if (leafHtSum == gct->getJetLeafCards().at(leaf)->getAllOutputHt().at(bx).value()) {
	htMinusVl += leafHtSum;
      } else { cout << "Ht sum check leaf " << leaf << endl; testPass = false; }
    }

    // Plus Wheel
    for ( ; leaf<6; leaf++) {
      unsigned leafHtSum = 0;
      for (int jf=0; jf<3; jf++) {
	assert (pJet != plusWheelJetData.end());
	leafHtSum += pJet->htSum;
	pJet++;
      }
      if (leafHtSum == gct->getJetLeafCards().at(leaf)->getAllOutputHt().at(bx).value()) {
	htPlusVal += leafHtSum;
      } else { cout << "Ht sum check leaf " << leaf << endl; testPass = false; }
    }

    unsigned htTotal = htMinusVl + htPlusVal;

    bool htMinusOvrFlow = (htMinusVl>=4096);
    bool htPlusOverFlow = (htPlusVal>=4096);
    bool htTotalOvrFlow = (htTotal>=4096) || htMinusOvrFlow  || htPlusOverFlow;
    //
    // Check the input to the final GlobalEnergyAlgos is as expected
    //--------------------------------------------------------------------------------------
    //
    if (!myGlobalEnergy->getInputHtVlMinusWheel().at(bx).overFlow() && !htMinusOvrFlow &&
	(myGlobalEnergy->getInputHtVlMinusWheel().at(bx).value()!=htMinusVl)) { cout << "ht Minus " << htMinusVl <<endl; testPass = false; }
    if (!myGlobalEnergy->getInputHtValPlusWheel().at(bx).overFlow() && !htPlusOverFlow &&
	(myGlobalEnergy->getInputHtValPlusWheel().at(bx).value()!=htPlusVal)) { cout << "ht Plus " << htPlusVal <<endl; testPass = false; }

    // Check the output value
    if (!myGlobalEnergy->getEtHadColl().at(bx).overFlow() && !htTotalOvrFlow &&
	(myGlobalEnergy->getEtHadColl().at(bx).value() != htTotal)) {cout << "Algo etHad" << endl; testPass = false;}
    // end of loop over bunch crossings
  }
  return testPass;
}

//=================================================================================================================
//
/// Check the jet counting algorithms
bool gctTestHtAndJetCounts::checkJetCounts(const L1GlobalCaloTrigger* gct) const
{
  bool testPass = true;
  L1GctGlobalEnergyAlgos* myGlobalEnergy = gct->getEnergyFinalStage();
  const L1GctJetEtCalibrationLut* myLut = gct->getJetEtCalibLut();

  int bx = m_bxStart;
  for (int ibx=0; ibx<m_numOfBx; ibx++) {
    //
    // Emulate the jet counting
    //--------------------------------------------------------------------------------------
    //
    vector<unsigned> JcMinusWheel(L1GctWheelJetFpga::N_JET_COUNTERS);
    vector<unsigned> JcPlusWheel (L1GctWheelJetFpga::N_JET_COUNTERS);
    vector<unsigned> JcResult(L1GctWheelJetFpga::N_JET_COUNTERS);

    for (unsigned jcnum = 0 ; jcnum<L1GctWheelJetFpga::N_JET_COUNTERS ; jcnum++) {

      L1GctJetCounterSetup::cutsListForJetCounter
	cutListPos=gct->getWheelJetFpgas().at(0)->getJetCounter(jcnum)->getJetCounterLut()->cutList();
      L1GctJetCounterSetup::cutsListForJetCounter
	cutListNeg=gct->getWheelJetFpgas().at(1)->getJetCounter(jcnum)->getJetCounterLut()->cutList();

      unsigned count0 = countJetsInCut(minusWheelJetDta, cutListNeg, myLut, bx) ;
      JcMinusWheel.at(jcnum) = count0;
      unsigned count1 = countJetsInCut(plusWheelJetData, cutListPos, myLut, bx) ;
      JcPlusWheel.at(jcnum) = count1;
      JcResult.at(jcnum) = ( (count0<7) && (count1<7) ? (count0 + count1) : 31 ) ;
    }

    // Check the inputs from the two wheels
    static const unsigned int N_JET_COUNTERS = L1GctWheelJetFpga::N_JET_COUNTERS;
    for (unsigned int i=0; i<N_JET_COUNTERS; i++) {
      if ((myGlobalEnergy->getInputJcVlMinusWheel().at(N_JET_COUNTERS*ibx + i).value()!=JcMinusWheel.at(i)) &&
	  (myGlobalEnergy->getInputJcVlMinusWheel().at(N_JET_COUNTERS*ibx + i).overFlow() ^ (JcMinusWheel.at(i)==7))) {
	cout << "jc Minus " << i << " value " << JcMinusWheel.at(i) <<endl;
	testPass = false;
      }
      if ((myGlobalEnergy->getInputJcValPlusWheel().at(N_JET_COUNTERS*ibx + i).value()!=JcPlusWheel.at(i)) ||
	  (myGlobalEnergy->getInputJcValPlusWheel().at(N_JET_COUNTERS*ibx + i).overFlow() ^ (JcPlusWheel.at(i)==7))) {
	cout << "jc Plus " << i << " value " << JcPlusWheel.at(i) <<endl;
	testPass = false;
      }
    }
    // Check the outputs
    for (unsigned int j=0 ; j<L1GctWheelJetFpga::N_JET_COUNTERS ; j++) {
      if ((myGlobalEnergy->getJetCount(j, bx).value() != JcResult.at(j)) || 
	  (myGlobalEnergy->getJetCount(j, bx).overFlow() ^ (JcResult.at(j)==31))) { 
	cout << "Algo jCount " << j << endl;
	cout << "Expected " << JcResult.at(j) << " found " << myGlobalEnergy->getJetCount(j, bx) << endl;
	cout << "PlusWheel " << myGlobalEnergy->getInputJcValPlusWheel().at(12*ibx + j) << endl;
	//      cout << *myGlobalEnergy->getPlusWheelJetFpga() << endl;
	cout << "MinusWheel " << myGlobalEnergy->getInputJcVlMinusWheel().at(12*ibx + j) << endl;
	//      cout << *myGlobalEnergy->getMinusWheelJetFpga() << endl;
	testPass = false;
      }
    }
    // end of loop over bunch crossings
    bx++;
  }
  return testPass;
}

//=================================================================================================================
//
// PRIVATE MEMBER FUNCTIONS
//
gctTestHtAndJetCounts::rawJetData gctTestHtAndJetCounts::rawJetFinderOutput(const L1GctJetFinderBase* jf, const int bx) const
{
  RawJetsVector jetsFromJf = jf->getRawJets();
  RawJetsVector jetList;
  unsigned sumHt = 0;
  for (RawJetsVector::const_iterator jet=jetsFromJf.begin(); jet!=jetsFromJf.end(); jet++) {
    if (jet->bx()==bx && !jet->isNullJet()) {
//       cout << "found a jet " << jet->rawsum()
// 	   << " eta " << jet->globalEta()
// 	   << " phi " << jet->globalPhi()
// 	   << " bx " << jet->bx() << endl;
      jetList.push_back(*jet);
      sumHt += jet->calibratedEt(jf->getJetEtCalLut());
    }
  }
  rawJetData result(jetList, sumHt);
  return result;
}

//
// Function definition for jet count checking
//=========================================================================
// Does what it says ...
unsigned gctTestHtAndJetCounts::countJetsInCut(const vector<rawJetData>& jetList,
                                               const L1GctJetCounterSetup::cutsListForJetCounter& cutList,
                                               const L1GctJetEtCalibrationLut* lut,
					       const int bx) const
{
  unsigned count = 0;
  L1GctJetCount<3> dummy;
  const unsigned MAX_VALUE = (1<<(dummy.size()))-1;

  // Loop over all jets in all jetFinders in a wheel
  for (vector<rawJetData>::const_iterator jList=jetList.begin(); jList<jetList.end(); jList++) {
    for (RawJetsVector::const_iterator jet=jList->jets.begin(); jet<jList->jets.end(); jet++) {
      if (jet->bx()==bx && !jet->jetCand(lut).empty()) {
        bool jetPassesCut = true;
        for (L1GctJetCounterSetup::cutsListForJetCounter::const_iterator cut=cutList.begin(); cut<cutList.end(); cut++) {
          switch (cut->cutType) {
            case L1GctJetCounterSetup::minRank :
              jetPassesCut &= (jet->jetCand(lut).rank() >= cut->cutValue1);
              break;

            case L1GctJetCounterSetup::maxRank:
              jetPassesCut &= (jet->jetCand(lut).rank() <= cut->cutValue1);
              break;

            case L1GctJetCounterSetup::centralEta:
              jetPassesCut &= (jet->rctEta() <= cut->cutValue1);
              break;

            case L1GctJetCounterSetup::forwardEta:
              jetPassesCut &= (jet->rctEta() >= cut->cutValue1);
              break;

            case L1GctJetCounterSetup::phiWindow:
              jetPassesCut &= (cut->cutValue2 > cut->cutValue1 ?
                ((jet->globalPhi() >= cut->cutValue1) && (jet->globalPhi() <= cut->cutValue2)) :
                ((jet->globalPhi() >= cut->cutValue1) || (jet->globalPhi() <= cut->cutValue2)));
              break;

            case L1GctJetCounterSetup::nullCutType:
              jetPassesCut &= false;
              break;

          }
        }
        if (jetPassesCut && (count<MAX_VALUE)) { count++; }
      }
    }
  }
  return count;
}



