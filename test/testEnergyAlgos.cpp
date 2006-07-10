/*! \file testEnergyAlgos.cpp
 * \brief Procedural skeleton unit-test code for the L1Gct(...)EnergyAlgos classes.
 *
 * Update comment when we have decided what the code does!
 *  This is skeleton code that reads in data from a text file to feed into
 *  the setInputRegions() method, runs the process() method, and then
 *  checks the data from the outputting methods getInputRegions() and
 *  getJets() against known results also stored in the text file.
 *
 * \author Greg Heath
 * \date March 2006
 */


//NOTE all these includes need to be sorted with a proper build include path,
//rather than a relative path, in order to comply with CMS style.

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"  //The class to be tested
#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <exception> //for exception handling
#include <stdexcept> //for std::runtime_error()
using namespace std;

//Typedefs for the vector templates used
        struct etmiss_vec { unsigned mag; unsigned phi;};
// typedef vector<L1CaloRegion> RegionsVector;
typedef vector<L1GctJet> JetsVector;

//  FUNCTION PROTOTYPES
/// Runs the test, and returns a string with the test result message in.
string classTest();
//
/// Generate an event, load it into the gct, and return some sums for checking
void loadNextEvent(L1GlobalCaloTrigger* &gct, const bool simpleEvent,
	       vector<unsigned> &etStripSums, bool &inMinusOvrFlow, bool &inPlusOverFlow);
/// Check the energy sums algorithms
bool checkEnergySums(const L1GlobalCaloTrigger* gct,
		     const vector<unsigned> &etStripSums,
		     const bool inMinusOvrFlow, const bool inPlusOverFlow);

/// Check the Ht summing algorithms
bool checkHtSums(const L1GlobalCaloTrigger* gct,
		 JetsVector &MinusWheelJets, JetsVector &PlusWheelJets);
/// Check the jet counting algorithms
bool checkJetCounts(const L1GlobalCaloTrigger* gct,
		    const JetsVector MinusWheelJets, const JetsVector PlusWheelJets);

/// Entrypoint of unit test code + error handling
int main(int argc, char **argv)
{
    cout << "\n*********************************" << endl;
    cout << "Hello from testEnergyAlgos       " << endl;
    cout << "*********************************" << endl;

    srand(290306);
    try
    {
        cout << "\n" << classTest() << endl;
    }
    catch(const exception &e)
    {
        cerr << "\nError! " << e.what() << endl;
    }
    catch(...)
    {
        cerr << "\nError! An unknown exception has occurred!" << endl;
    }
    return 0;   
}

// Runs the test, and returns a string with the test result message in.
string classTest()
{
  bool testPass = true;       //Test passing flag.
    
  const int maxTests=10000;
  const int initialTests=100;
 
  // Initialise the gct
  L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(false, L1GctJetLeafCard::hardwareJetFinder);
//   L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(false, L1GctJetLeafCard::tdrJetFinder);

  for (int t=0; t<maxTests; t++)
    {
      if (t <initialTests) { cout << "Test " << t << endl; }
      if (t==initialTests) { cout << "Test class has passed initial simple tests!" << endl; }
      // Initialise the gct
      gct->reset();

      // Generate an event
      vector <unsigned>etStripSums(36);
      bool inMinusOvrFlow, inPlusOverFlow;
      loadNextEvent(gct, (t<initialTests), etStripSums, inMinusOvrFlow, inPlusOverFlow);

      //
      // Run the processing
      //--------------------------------------------------------------------------------------
      //
   
      gct->process();  //Run algorithm

      //
      // Check the energy sum logic
      //--------------------------------------------------------------------------------------
      //

      testPass &= checkEnergySums(gct, etStripSums,inMinusOvrFlow, inPlusOverFlow);
      if(testPass == false)
	{
	  // Print failed events for debug purposes
	  cout << *gct->getEnergyFinalStage() << endl;
	  cout << "Failed at test number " << t <<endl;
	  return "Test class has failed check of energy summing logic!";
	}

      //
      // Check the Ht summing
      //--------------------------------------------------------------------------------------
      //
    
      // This check fills lists of the jets found in the event,
      // which are then used by the jet count checking
      JetsVector MinusWheelJets;
      JetsVector PlusWheelJets;

      testPass &= checkHtSums(gct, MinusWheelJets, PlusWheelJets); 
      if(testPass == false)
	{
	  // Print failed events for debug purposes
	  cout << *gct->getEnergyFinalStage() << endl;
	  cout << "Failed at test number " << t <<endl;
	  return "Test class has failed check of Ht summing logic!";
	}
   
      //
      // Check the jet counts
      //--------------------------------------------------------------------------------------
      //
    
      testPass &= checkJetCounts(gct, MinusWheelJets, PlusWheelJets); 
      if(testPass == false)
	{
	  // Print failed events for debug purposes
	  cout << *gct->getEnergyFinalStage() << endl;
	  cout << "Failed at test number " << t <<endl;
	  return "Test class has failed check of jet count logic!";
	}
    }
  // Test the print routine at the end
  cout << *gct->getEnergyFinalStage() << endl;
  return "Test class has passed!";
}

//
//=========================================================================
// Here's the event generation
//=========================================================================
//
// FUNCTION PROTOTYPES FOR EVENT GENERATION
/// Generates test data for missing Et as 2-vector (magnitude, direction)
etmiss_vec randomMissingEtVector();
/// Generates test data consisting of energies to be added together with their sum
void generateTestData(vector<unsigned> &energies, int size, unsigned max);
//=========================================================================
/// Generate an event, load it into the gct, and return some sums for checking
void loadNextEvent(L1GlobalCaloTrigger* &gct, const bool simpleEvent,
	       vector<unsigned> &etStripSums, bool &inMinusOvrFlow, bool &inPlusOverFlow)
{
  for (int i=0; i<36; i++) {
    etStripSums.at(i)=0;
  }
  inMinusOvrFlow = false;
  inPlusOverFlow = false;

  // For initial tests just try things out with one region input
  // Then test with summing multiple regions. Choose one value
  // of energy and phi for each eta to avoid trying to set the
  // same region several times.
  for (unsigned i=0; i<(simpleEvent ? 1 : L1CaloRegionDetId::N_ETA); i++) {
    etmiss_vec etVector=randomMissingEtVector();
//     cout << "Region et " << etVector.mag << " phi " << etVector.phi << endl;
    // Set a single region input
    unsigned etaRegion = i;
    unsigned phiRegion = etVector.phi/4;

    gct->setRegion(etVector.mag, etaRegion, phiRegion);
        
    // Here we fill the expected values. Et values restricted to
    // eight bits in HF and ten bits in the rest of the system.
    if (etaRegion<(L1CaloRegionDetId::N_ETA)/2) {
      if (etaRegion<4) {
	etStripSums.at(phiRegion) += (etVector.mag & 0xff);
	inMinusOvrFlow |= (etVector.mag>=0x100);
      } else {
	etStripSums.at(phiRegion) += (etVector.mag & 0x3ff);
	inMinusOvrFlow |= (etVector.mag>=0x400);
      }
    } else {
      if (etaRegion>=18) {
	etStripSums.at(phiRegion+L1CaloRegionDetId::N_PHI) += (etVector.mag & 0xff);
	inPlusOverFlow |= (etVector.mag>=0x100);
      } else {
	etStripSums.at(phiRegion+L1CaloRegionDetId::N_PHI) += (etVector.mag & 0x3ff);
	inPlusOverFlow |= (etVector.mag>=0x400);
      }
    }
  }
}

//
// Function definitions for event generation
//=========================================================================
// Generates 2-d missing Et vector
etmiss_vec randomMissingEtVector()
{
  // This produces random variables distributed as a 2-d Gaussian
  // with a standard deviation of sigma for each component,
  // and magnitude ranging up to 5*sigma.
  //
  // With sigma set to 400 we will always be in the range
  // of 12-bit signed integers (-2048 to 2047).
  // With sigma set to 200 we are always in the range
  // of 10-bit input region Et values.
  const float sigma=200.;

  // rmax controls the magnitude range
  // Chosen as a power of two conveniently close to
  // exp(5*5/2) to give a 5*sigma range.
  const unsigned rmax=262144;

  vector<unsigned> components(2);
  float p,r,s;
  unsigned Emag, Ephi;

  const float nbins = 18.;

  // Generate a pair of uniform pseudo-random integers
  generateTestData(components, (int) 2, rmax);

  // Exclude the value zero for the first random integer
  // (Alternatively, return an overflow bit)
  while (components[0]==0) {generateTestData(components, (int) 2, rmax);}

  // Convert to the 2-d Gaussian
  r = float(rmax);
  s = r/float(components[0]);
  p = float(components[1])/r;
  // Force phi value into the centre of a bin
  Emag = int(sigma*sqrt(2.*log(s)));
  Ephi = int(nbins*p);
  // Copy to the output
  etmiss_vec Et;
  Et.mag = Emag;
  Et.phi = (4*Ephi);
  return Et;

}

/// Generates test data consisting of a vector of energies
/// uniformly distributed between zero and max
void generateTestData(vector<unsigned> &energies, int size, unsigned max)
{
  int i;
  int r,e;
  float p,q,s,t;

  p = float(max);
  q = float(RAND_MAX);
  for (i=0; i<size; i++) {
    r = rand();
    s = float(r);
    t = s*p/q;
    e = int(t);

    energies[i] = e;
  }
}


//
//=========================================================================
// Here's the procedure for checking the energy sums (Ex, Ey, total and missing Et)
//=========================================================================
//
// FUNCTION PROTOTYPES FOR ENERGY SUM CHECKING
/// Integer calculation of Ex or Ey from magnitude for a given phi bin
int etComponent(const unsigned Emag, const unsigned fact);
/// Calculate et vector from ex and ey, using floating arithmetic and conversion back to integer
etmiss_vec trueMissingEt(const int ex, const int ey);
//=========================================================================
/// Check the energy sums algorithms
bool checkEnergySums(const L1GlobalCaloTrigger* gct,
		     const vector<unsigned> &etStripSums,
		     const bool inMinusOvrFlow, const bool inPlusOverFlow)
{
  bool testPass=true;
  L1GctGlobalEnergyAlgos* myGlobalEnergy = gct->getEnergyFinalStage();

  int exMinusVl = 0;
  int eyMinusVl = 0;
  unsigned etMinusVl = 0;

  int strip = 0;
  // Find the expected sums for the Minus end
  for ( ; strip<18; strip++) {
    unsigned et = etStripSums.at(strip);

    unsigned rctStrip = (40-strip)%18 + 18*(strip/18);
    L1GctJetLeafCard* jlc = gct->getJetLeafCards().at(rctStrip/6);
    L1GctJetFinderBase* jf = ((rctStrip%6)/2==0) ? jlc->getJetFinderA() :
      ( ((rctStrip%6)/2==1) ? jlc->getJetFinderB() : jlc->getJetFinderC() );
    L1GctScalarEtVal gctEt = ((rctStrip%2==0) ? jf->getEtStrip0() : jf->getEtStrip1() );

    int ex = etComponent(et, ((2*strip+9)%36) );
    int ey = etComponent(et, (( 2*strip )%36) );

    exMinusVl += ex;
    eyMinusVl += ey;
    etMinusVl += et; 
  }
  bool exMinusOvrFlow = (exMinusVl<-2048) || (exMinusVl>=2048) || inMinusOvrFlow;
  bool eyMinusOvrFlow = (eyMinusVl<-2048) || (eyMinusVl>=2048) || inMinusOvrFlow;
  bool etMinusOvrFlow = (etMinusVl>=4096) || inMinusOvrFlow;

  int exPlusVal = 0;
  int eyPlusVal = 0;
  unsigned etPlusVal = 0;

  // Find the expected sums for the Plus end
  for ( ; strip<36; strip++) {
    unsigned et = etStripSums.at(strip);

    unsigned rctStrip = (40-strip)%18 + 18*(strip/18);
    L1GctJetLeafCard* jlc = gct->getJetLeafCards().at(rctStrip/6);
    L1GctJetFinderBase* jf = ((rctStrip%6)/2==0) ? jlc->getJetFinderA() :
      ( ((rctStrip%6)/2==1) ? jlc->getJetFinderB() : jlc->getJetFinderC() );
    L1GctScalarEtVal gctEt = ((rctStrip%2==0) ? jf->getEtStrip0() : jf->getEtStrip1() );

    int ex = etComponent(et, ((2*strip+9)%36) );
    int ey = etComponent(et, (( 2*strip )%36) );

    exPlusVal += ex;
    eyPlusVal += ey;
    etPlusVal += et; 
  }
  bool exPlusOverFlow = (exPlusVal<-2048) || (exPlusVal>=2048) || inPlusOverFlow;
  bool eyPlusOverFlow = (eyPlusVal<-2048) || (eyPlusVal>=2048) || inPlusOverFlow;
  bool etPlusOverFlow = (etPlusVal>=4096) || inPlusOverFlow;

  int exTotal = exMinusVl + exPlusVal;
  int eyTotal = eyMinusVl + eyPlusVal;
  unsigned etTotal = etMinusVl + etPlusVal;

  etmiss_vec etResult = trueMissingEt(exTotal, eyTotal);

  bool exTotalOvrFlow = (exTotal<-2048) || (exTotal>=2048) || exMinusOvrFlow || exPlusOverFlow;
  bool eyTotalOvrFlow = (eyTotal<-2048) || (eyTotal>=2048) || eyMinusOvrFlow || eyPlusOverFlow;
  bool etTotalOvrFlow = (etTotal>=4096) || etMinusOvrFlow  || etPlusOverFlow;

  bool etMissOverFlow = exTotalOvrFlow || eyTotalOvrFlow;

  //
  // Check the input to the final GlobalEnergyAlgos is as expected
  //--------------------------------------------------------------------------------------
  //
  if (!myGlobalEnergy->getInputExVlMinusWheel().overFlow() && !exMinusOvrFlow &&
      (myGlobalEnergy->getInputExVlMinusWheel().value()!=exMinusVl)) { cout << "ex Minus " << exMinusVl <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputExValPlusWheel().overFlow() && !exPlusOverFlow &&
      (myGlobalEnergy->getInputExValPlusWheel().value()!=exPlusVal)) { cout << "ex Plus " << exPlusVal <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputEyVlMinusWheel().overFlow() && !eyMinusOvrFlow &&
      (myGlobalEnergy->getInputEyVlMinusWheel().value()!=eyMinusVl)) { cout << "ey Minus " << eyMinusVl <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputEyValPlusWheel().overFlow() && !eyPlusOverFlow &&
      (myGlobalEnergy->getInputEyValPlusWheel().value()!=eyPlusVal)) { cout << "ey Plus " << eyPlusVal <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputEtVlMinusWheel().overFlow() && !etMinusOvrFlow &&
      (myGlobalEnergy->getInputEtVlMinusWheel().value()!=etMinusVl)) { cout << "et Minus " << etMinusVl <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputEtValPlusWheel().overFlow() && !etPlusOverFlow &&
      (myGlobalEnergy->getInputEtValPlusWheel().value()!=etPlusVal)) { cout << "et Plus " << etPlusVal <<endl; testPass = false; }
 
  //
  // Now check the processing in the final stage GlobalEnergyAlgos
  //--------------------------------------------------------------------------------------
  //
  // Check the missing Et calculation. Allow some margin for the
  // integer calculation of missing Et.
  if (!etMissOverFlow && !myGlobalEnergy->getEtMiss().overFlow()) {
    unsigned etDiff, phDiff;
    unsigned etMargin, phMargin;

    etDiff = (unsigned) abs((long int) etResult.mag - (long int) myGlobalEnergy->getEtMiss().value());
    phDiff = (unsigned) abs((long int) etResult.phi - (long int) myGlobalEnergy->getEtMissPhi().value());
    if (phDiff>60) {phDiff=72-phDiff;}
    //
    etMargin = max((etResult.mag/100), (unsigned) 1) + 2;
    if (etResult.mag==0) { phMargin = 72; } else { phMargin = (30/etResult.mag) + 1; }
    if ((etDiff > etMargin) || (phDiff > phMargin)) {cout << "Algo etMiss diff "
							  << etDiff << " phi diff " << phDiff << endl; testPass = false;}
  }
  // Check the total Et calculation
  if (!myGlobalEnergy->getEtSum().overFlow() && !etTotalOvrFlow &&
      (myGlobalEnergy->getEtSum().value() != etTotal)) {cout << "Algo etSum" << endl; testPass = false;}
  return testPass;
}

//
// Function definitions for energy sum checking
//=========================================================================
int etComponent(const unsigned Emag, const unsigned fact) {
  // Copy the Ex, Ey conversion from the hardware emulation
  const unsigned sinFact[10] = {0, 44, 87, 128, 164, 196, 221, 240, 252, 256};
  unsigned myFact;
  bool negativeResult;
  int result;
  switch (fact/9) {
  case 0:
    myFact = sinFact[fact];
    negativeResult = false;
    break;
  case 1:
    myFact = sinFact[(18-fact)];
    negativeResult = false;
    break;
  case 2:
    myFact = sinFact[(fact-18)];
    negativeResult = true;
    break;
  case 3:
    myFact = sinFact[(36-fact)];
    negativeResult = true;
    break;
  default:
    cout << "Invalid factor " << fact << endl;
    return 0;
  }
  result = static_cast<int>(Emag*myFact);
  // Divide by 256 using bit-shift; but emulate
  // twos-complement arithmetic for negative numbers
  if ( negativeResult ) {
    result = (1<<24)-result;
    result = result>>8;
    result = result-(1<<16);
  } else { result = result>>8; }
  return result;
}

/// Calculate the expected missing Et vector for a given
/// ex and ey sum, for comparison with the hardware
etmiss_vec trueMissingEt(const int ex, const int ey) {

  etmiss_vec result;

  double fx = static_cast<double>(ex);
  double fy = static_cast<double>(ey);
  double fmag = sqrt(fx*fx + fy*fy);
  double fphi = 36.*atan2(fy, fx)/3.1415927;

  result.mag = static_cast<int>(fmag);
  if (fphi>=0) {
    result.phi = static_cast<int>(fphi);
  } else {
    result.phi = static_cast<int>(fphi+72.);
  }

  return result;

}

//
//=========================================================================
// Here's the procedure for checking the jet-energy sum (Ht calculation)
//=========================================================================
//
// FUNCTION PROTOTYPES FOR HT SUM CHECKING
/// Checks the Ht calculation in a leaf card and returns the Ht sum
bool checkHt(L1GctJetLeafCard* jlc, JetsVector &jetList, unsigned &leafHt);
/// Works out the Ht for a jet from the raw input regions
unsigned jetHtSum(L1GctJetFinderBase* jf, int jn);
//=========================================================================
/// Check the Ht summing algorithms
bool checkHtSums(const L1GlobalCaloTrigger* gct, JetsVector &MinusWheelJets, JetsVector &PlusWheelJets)
{
  bool testPass = true;
  L1GctGlobalEnergyAlgos* myGlobalEnergy = gct->getEnergyFinalStage();
  //
  // Check the Ht calculation (starting from the found jets)
  //--------------------------------------------------------------------------------------
  //
  // Minus Wheel
  unsigned htMinusVl = 0;
  MinusWheelJets.clear();
  for (int leaf=0; leaf<3; leaf++) {
    unsigned ht;
    if (checkHt(gct->getJetLeafCards().at(leaf), MinusWheelJets, ht)) {
      htMinusVl += ht;
    } else { cout << "Ht sum check leaf " << leaf << endl; testPass = false; }
  }

  // Plus Wheel
  unsigned htPlusVal = 0;
  PlusWheelJets.clear();
  for (int leaf=3; leaf<6; leaf++) {
    unsigned ht = 0;
    if (checkHt(gct->getJetLeafCards().at(leaf), PlusWheelJets, ht)) {
      htPlusVal += ht;
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
  if (!myGlobalEnergy->getInputHtVlMinusWheel().overFlow() && !htMinusOvrFlow &&
      (myGlobalEnergy->getInputHtVlMinusWheel().value()!=htMinusVl)) { cout << "ht Minus " << htMinusVl <<endl; testPass = false; }
  if (!myGlobalEnergy->getInputHtValPlusWheel().overFlow() && !htPlusOverFlow &&
      (myGlobalEnergy->getInputHtValPlusWheel().value()!=htPlusVal)) { cout << "ht Plus " << htPlusVal <<endl; testPass = false; }

  // Check the output value
  if (!myGlobalEnergy->getEtHad().overFlow() && !htTotalOvrFlow &&
      (myGlobalEnergy->getEtHad().value() != htTotal)) {cout << "Algo etHad" << endl; testPass = false;}
  return testPass;
}

/// Check the initial Ht calculation
bool checkHt(L1GctJetLeafCard* jlc, JetsVector &jetList, unsigned &leafHt) {

  unsigned result=jlc->getOutputHt().value();
  unsigned sumHt=0;

  for (unsigned i=0; i<L1GctJetFinderBase::MAX_JETS_OUT; i++) {
    if (!jlc->getOutputJetsA().at(i).isNullJet()) {
//       cout << "found a jet " << jlc->getOutputJetsA().at(i).rank()
// 	   << " eta " << jlc->getOutputJetsA().at(i).globalEta()
// 	   << " phi " << jlc->getOutputJetsA().at(i).globalPhi() << endl;
      jetList.push_back(jlc->getOutputJetsA().at(i));
      sumHt += jetHtSum(jlc->getJetFinderA(), i);
    }
    if (!jlc->getOutputJetsB().at(i).isNullJet()) {
//       cout << "found a jet " << jlc->getOutputJetsB().at(i).rank()
// 	   << " eta " << jlc->getOutputJetsB().at(i).globalEta()
// 	   << " phi " << jlc->getOutputJetsB().at(i).globalPhi() << endl;
      jetList.push_back(jlc->getOutputJetsB().at(i));
      sumHt += jetHtSum(jlc->getJetFinderB(), i);
    }
    if (!jlc->getOutputJetsC().at(i).isNullJet()) {
//       cout << "found a jet " << jlc->getOutputJetsC().at(i).rank()
// 	   << " eta " << jlc->getOutputJetsC().at(i).globalEta()
// 	   << " phi " << jlc->getOutputJetsC().at(i).globalPhi() << endl;
      jetList.push_back(jlc->getOutputJetsC().at(i));
      sumHt += jetHtSum(jlc->getJetFinderC(), i);
    }
  }

  leafHt = result;
  return (sumHt==result);
}

//
// Function definitions for Ht sum checking
//=========================================================================
/// Work out the Ht for a jet.
unsigned jetHtSum(L1GctJetFinderBase* jf, int jn) {

  //
  return static_cast<unsigned>(jf->getJets().at(jn).rankForHt());
  //

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  // // The following code (commented) goes back to the original
  // // inputRegions and works out the expected jet rank.
  // // It only works for the tdr jetFinder though.
//   vector<L1CaloRegion>inputRegions = jf->getInputRegions();
//   const unsigned COL_OFFSET = ((L1CaloRegionDetId::N_ETA)/2)+1;

//   // Check the input array size
//   if (inputRegions.size()!=(COL_OFFSET*4)) {
//     cout << "Invalid size for jet finder input " << inputRegions.size() << " expecting " << (COL_OFFSET*4) << endl;
//     return 0;
//   }

//   // Find the eta and phi for this jet, and check
//   unsigned eta = static_cast<unsigned>(jf->getJets().at(jn).rctEta());
//   unsigned phi = static_cast<unsigned>(jf->getJets().at(jn).rctPhi());

//   if (phi>1 || eta>=(COL_OFFSET-1)) {
//     cout << "Invalid eta, phi for jet: " << eta << ", " << phi << endl;
//     return 0;
//   }

//   // Sum the et values for the nine regions centred on this eta and phi
//   unsigned rawSum = 0;

//   for (unsigned col=phi; col<=phi+2; col++) {
//     rawSum += inputRegions.at(col*COL_OFFSET+eta).et();
//     rawSum += inputRegions.at(col*COL_OFFSET+eta+1).et();
//     if ((eta+2)<COL_OFFSET) {
//       rawSum += inputRegions.at(col*COL_OFFSET+eta+2).et();
//     }
//   }
//   // Convert to Ht and return
//   return static_cast<unsigned>(jf->getJetEtCalLut()->convertToTenBitRank(static_cast<uint16_t>(rawSum),
// 									 static_cast<uint16_t>(eta)));
  // // end of commented code
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
}

//
//=========================================================================
// Here's the procedure for checking the jet counting
//=========================================================================
//
// FUNCTION PROTOTYPE FOR JET COUNTING
/// Counts jets in cuts
unsigned countJetsInCut(const JetsVector &jetList, const unsigned jcnum, const unsigned Wheel);
//=========================================================================
/// Check the jet counting algorithms
bool checkJetCounts(const L1GlobalCaloTrigger* gct,
		    const JetsVector MinusWheelJets, const JetsVector PlusWheelJets)
{
  bool testPass = true;
  L1GctGlobalEnergyAlgos* myGlobalEnergy = gct->getEnergyFinalStage();
  //
  // Emulate the jet counting
  //--------------------------------------------------------------------------------------
  //
  vector<unsigned> JcMinusWheel(L1GctWheelJetFpga::N_JET_COUNTERS);
  vector<unsigned> JcPlusWheel (L1GctWheelJetFpga::N_JET_COUNTERS);
  vector<unsigned> JcResult(L1GctWheelJetFpga::N_JET_COUNTERS);

  for (unsigned jcnum = 0 ; jcnum<L1GctWheelJetFpga::N_JET_COUNTERS ; jcnum++) {
    unsigned count0 = countJetsInCut(MinusWheelJets, jcnum, 0) ;
    JcMinusWheel.at(jcnum) = count0;
    unsigned count1 = countJetsInCut(PlusWheelJets, jcnum, 1) ;
    JcPlusWheel.at(jcnum) = count1;
    JcResult.at(jcnum) = ( (count0<7) && (count1<7) ? (count0 + count1) : 31 ) ;
  }

  // Check the inputs from the two wheels
  for (unsigned int i=0; i<L1GctWheelJetFpga::N_JET_COUNTERS; i++) {
    if ((myGlobalEnergy->getInputJcVlMinusWheel(i).value()!=JcMinusWheel.at(i)) &&
	(myGlobalEnergy->getInputJcVlMinusWheel(i).overFlow() ^ (JcMinusWheel.at(i)==7))) {
      cout << "jc Minus " << i << " value " << JcMinusWheel.at(i) <<endl;
      testPass = false;
    }
    if ((myGlobalEnergy->getInputJcValPlusWheel(i).value()!=JcPlusWheel.at(i)) ||
	(myGlobalEnergy->getInputJcValPlusWheel(i).overFlow() ^ (JcPlusWheel.at(i)==7))) {
      cout << "jc Plus " << i << " value " << JcPlusWheel.at(i) <<endl;
      testPass = false;
    }
  }

  // Check the outputs
  for (unsigned int j=0 ; j<L1GctWheelJetFpga::N_JET_COUNTERS ; j++) {
    if ((myGlobalEnergy->getJetCount(j).value() != JcResult.at(j)) || 
	(myGlobalEnergy->getJetCount(j).overFlow() ^ (JcResult.at(j)==31))) { 
      cout << "Algo jCount " << j << endl;
      cout << "Expected " << JcResult.at(j) << " found " << myGlobalEnergy->getJetCount(j) << endl;
      cout << "PlusWheel " << myGlobalEnergy->getInputJcValPlusWheel(j) << endl;
      cout << *myGlobalEnergy->getPlusWheelJetFpga() << endl;
      cout << "MinusWheel " << myGlobalEnergy->getInputJcVlMinusWheel(j) << endl;
      cout << *myGlobalEnergy->getMinusWheelJetFpga() << endl;
      testPass = false;
    }
  }
  return testPass;
}

//
// Function definition for jet count checking
//=========================================================================
// Does what it says ...
unsigned countJetsInCut(const JetsVector &jetList, const unsigned jcnum, const unsigned Wheel)
{
  unsigned count = 0;
  L1GctJcWheelType dummy;
  const unsigned MAX_VALUE = (1<<(dummy.size()))-1;

  for (unsigned i=0; i<jetList.size(); i++) {
    bool jetPassesCut = false;
    switch (jcnum) {
    case (0) :
      jetPassesCut = (jetList.at(i).rank() >= 5);
      break;

    case (1) :
      jetPassesCut = (jetList.at(i).globalEta() >=5) && (jetList.at(i).globalEta() <= 16);
      break;

    case (2) :
      jetPassesCut = (jetList.at(i).globalEta() < 5) && (Wheel == 0);
      break;

    case (3) :
      jetPassesCut = (jetList.at(i).globalEta() >16) && (Wheel == 1);
      break;

    default :
      break;

    }
    if (jetPassesCut && (count<MAX_VALUE)) { count++; }
  }
  return count;
}

// // Loads test input regions from a text file.
// void loadTestData(RegionsVector &regions, JetsVector &jets, const string &fileName)
// {
//     // File input stream
//     ifstream fin;
    
//     safeOpenFile(fin, fileName);  //open the file
    
//     unsigned long int tempEt = 0;
//     unsigned short int tempMip = 0;
//     unsigned short int tempQuiet = 0;
    
//     // Loads the input data
//     for(unsigned int i = 0; i < regions.size(); ++i)
//     {
//         //read in the data from the file
//         fin >> tempEt;
//         fin >> tempMip;
//         fin >> tempQuiet;
        
//         regions[i].setEt(tempEt);
//         if(tempMip == 0) { regions[i].setMip(false); } else { regions[i].setMip(true); }
//         if(tempQuiet == 0) { regions[i].setQuiet(false); } else { regions[i].setQuiet(true); }
//     }
    
//     // Do similar to load the 'known' output jets (that we currently don't know...)
    
//     // Close the file
//     fin.close();    
        
//     return;
// }
    
    
// // Function to safely open files of any name, using a referenced return ifstream
// void safeOpenFile(ifstream &fin, const string &name)
// {
//     //Opens the file
//     fin.open(name.c_str(), ios::in);

//     //Error message, and return false if it goes pair shaped
//     if(!fin.good())
//     {
//         throw std::runtime_error("Couldn't open the file " + name + "!");
//     }
//     return;
// }
