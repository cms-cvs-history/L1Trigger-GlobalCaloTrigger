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
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinder.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetLeafCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetEtCalibrationLut.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h"

//Custom headers needed for this test
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctMap.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctRegion.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetCand.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEtTypes.h"

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
// typedef vector<L1GctRegion> RegionsVector;
// typedef vector<L1GctJet> JetsVector;

//  FUNCTION PROTOTYPES
/// Runs the test, and returns a string with the test result message in.
string classTest();
//
/// Generates test data for missing Et both as (x,y) components and
/// 2-vector (magnitude, direction)
void generateMissingEtTestData(int &Ex, int &Ey, etmiss_vec &Et);
/// Integer calculation of Ex or Ey from magnitude for a given phi bin
int etComponent(const unsigned Emag, const unsigned fact);
/// Calculate et vector from ex and ey, using floating arithmetic and conversion back to integer
etmiss_vec trueMissingEt(const int ex, const int ey);
/// Generates test data consisting of energies to be added together with their sum
void generateTestData(vector<unsigned> &energies, int size, unsigned max, unsigned &sum);
/// Checks the Ht calculation in a leaf card and returns the Ht sum
bool checkHt(L1GctJetLeafCard* jlc, unsigned &leafHt);
/// Works out the Ht for a jet from the raw input regions
unsigned jetHtSum(L1GctJetFinder* jf, int jn);

// /// Loads test input regions and also the known results from a text file.
// void loadTestData(RegionsVector &regions, JetsVector &jets, const string &fileName);
// /// Function to safely open files of any name, using a referenced return ifstream
// void safeOpenFile(ifstream &fin, const string &name);

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
  L1GctGlobalEnergyAlgos* myGlobalEnergy;

  bool testPass = true;       //Test passing flag.
    
  const int maxValues=3;
  vector<unsigned> energyValues(maxValues);
  vector<unsigned> JcResult(12);

  vector <unsigned>etStripSums(36);

  int exPlusVal, exMinusVl;
  int eyPlusVal, eyMinusVl;
  unsigned etPlusVal, etMinusVl;
  unsigned htPlusVal, htMinusVl;
  int exValue, eyValue;
  int exTotal, eyTotal;
  unsigned etTotal, htTotal;
  etmiss_vec etVector, etResult;

  bool inPlusOverFlow, inMinusOvrFlow;
  bool exPlusOverFlow, eyPlusOverFlow, etPlusOverFlow, htPlusOverFlow;
  bool exMinusOvrFlow, eyMinusOvrFlow, etMinusOvrFlow, htMinusOvrFlow;
  bool exTotalOvrFlow, eyTotalOvrFlow, etTotalOvrFlow, htTotalOvrFlow;
  bool etMissOverFlow;

  unsigned etDiff, phDiff;
  unsigned etMargin, phMargin;

  const int maxTests=10000;
  const int initialTests=100;
 
  // For (eta,phi) mapping of input regions
  L1GctMap* map = L1GctMap::getMap();

  // Initialise the gct
  L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(false);
  myGlobalEnergy = gct->getEnergyFinalStage();

  for (int t=0; t<maxTests; t++)
    {
      if (t<initialTests) { cout << "Test " << t << endl; }
      // Initialise the gct
      gct->reset();

      for (int i=0; i<36; i++) { etStripSums[i]=0; }

      exPlusVal = 0;
      eyPlusVal = 0;
      etPlusVal = 0;
      htPlusVal = 0;
      exMinusVl = 0;
      eyMinusVl = 0;
      etMinusVl = 0;
      htMinusVl = 0;

      inPlusOverFlow = false;
      inMinusOvrFlow = false;
      htPlusOverFlow = false;
      htMinusOvrFlow = false;

      etResult.mag = 0;
      etResult.phi = 0;
      // For initial tests just try things out with one region input
      // Then test with summing multiple regions. Choose one value
      // of energy and phi for each eta to avoid trying to set the
      // same region several times.
      for (int i=0; i<(t<initialTests ? 1 : 22); i++) {
        generateMissingEtTestData(exValue, eyValue, etVector);
        //
        // Set a single region input
	unsigned etaRegion = i;
        unsigned phiRegion = etVector.phi/4;
        
        L1GctRegion temp(map->id(etaRegion,phiRegion), etVector.mag, false, false, false, false);
	gct->setRegion(temp);

        // Here we fill the expected values
	if (etaRegion<11) {
	  etStripSums[phiRegion] += etVector.mag;
	  inMinusOvrFlow |= (etVector.mag>=0x400);
	} else {
	  etStripSums[phiRegion+18] += etVector.mag;
	  inPlusOverFlow |= (etVector.mag>=0x400);
	}

      }

      // Find the expected sums
      int strip = 0;
      for (int leaf=0; leaf<3; leaf++) {
	for (int i=0; i<6; i++) {
	  unsigned et = etStripSums.at(strip);
	  int ex = etComponent(et, ((2*strip+10)%36) );
	  int ey = etComponent(et, ((2*strip+1 )%36) );
	  strip++;

	  exMinusVl += ex;
	  eyMinusVl += ey;
	  etMinusVl += et; 
	}
      }
      exMinusOvrFlow = (exMinusVl<-2048) || (exMinusVl>=2048) || inMinusOvrFlow;
      eyMinusOvrFlow = (eyMinusVl<-2048) || (eyMinusVl>=2048) || inMinusOvrFlow;
      etMinusOvrFlow = (etMinusVl>=4096) || inMinusOvrFlow;

      for (int leaf=3; leaf<6; leaf++) {
	for (int i=0; i<6; i++) {
	  unsigned et = etStripSums.at(strip);
	  int ex = etComponent(et, ((2*strip+10)%36) );
	  int ey = etComponent(et, ((2*strip+1 )%36) );
	  strip++;

	  exPlusVal += ex;
	  eyPlusVal += ey;
	  etPlusVal += et; 
	}
      }
      exPlusOverFlow = (exPlusVal<-2048) || (exPlusVal>=2048) || inPlusOverFlow;
      eyPlusOverFlow = (eyPlusVal<-2048) || (eyPlusVal>=2048) || inPlusOverFlow;
      etPlusOverFlow = (etPlusVal>=4096) || inPlusOverFlow;

      exTotal = exMinusVl + exPlusVal;
      eyTotal = eyMinusVl + eyPlusVal;
      etTotal = etMinusVl + etPlusVal;

      etResult = trueMissingEt(exTotal, eyTotal);

      exTotalOvrFlow = (exTotal<-2048) || (exTotal>=2048) || exMinusOvrFlow || exPlusOverFlow;
      eyTotalOvrFlow = (eyTotal<-2048) || (eyTotal>=2048) || eyMinusOvrFlow || eyPlusOverFlow;
      etTotalOvrFlow = (etTotal>=4096) || etMinusOvrFlow  || etPlusOverFlow;

      etMissOverFlow = exTotalOvrFlow || eyTotalOvrFlow;

      //
      // Run the processing
      //--------------------------------------------------------------------------------------
      //
   
      gct->process();  //Run algorithm

      //
      // Check the Ht calculation (starting from the found jets)
      //--------------------------------------------------------------------------------------
      //
      // Minus Wheel
      for (int leaf=0; leaf<3; leaf++) {
	unsigned ht;
	if (checkHt(gct->getJetLeafCards().at(leaf), ht)) {
	  htMinusVl += ht;
	} else { cout << "Ht sum check leaf " << leaf << endl; testPass = false; }
      }
      // Plus Wheel
      for (int leaf=3; leaf<6; leaf++) {
	unsigned ht = 0;
	if (checkHt(gct->getJetLeafCards().at(leaf), ht)) {
	  htPlusVal += ht;
	} else { cout << "Ht sum check leaf " << leaf << endl; testPass = false; }
      }
      htMinusOvrFlow = (htMinusVl>=4096);
      htPlusOverFlow = (htPlusVal>=4096);
      htTotal = htMinusVl + htPlusVal;
      htTotalOvrFlow = (htTotal>=4096) || htMinusOvrFlow  || htPlusOverFlow;
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
      if (!myGlobalEnergy->getInputHtVlMinusWheel().overFlow() && !htMinusOvrFlow &&
	  (myGlobalEnergy->getInputHtVlMinusWheel().value()!=htMinusVl)) { cout << "ht Minus " << htMinusVl <<endl; testPass = false; }
      if (!myGlobalEnergy->getInputHtValPlusWheel().overFlow() && !htPlusOverFlow &&
	  (myGlobalEnergy->getInputHtValPlusWheel().value()!=htPlusVal)) { cout << "ht Plus " << htPlusVal <<endl; testPass = false; }
   
      if(testPass == false)
	{
	  cout << *myGlobalEnergy << endl;
	  cout << "Failed at test number " << t <<endl;
	  return "Test class has failed initial data input/output comparison!";
	}
      //
      //--------------------------------------------------------------------------------------
      //
      // Check the missing Et calculation. Allow some margin for the
      // integer calculation of missing Et.
      if (!etMissOverFlow && !myGlobalEnergy->getEtMiss().overFlow()) {
        etDiff = (unsigned) abs((long int) etResult.mag - (long int) myGlobalEnergy->getEtMiss().value());
        phDiff = (unsigned) abs((long int) etResult.phi - (long int) myGlobalEnergy->getEtMissPhi().value());
        if (phDiff>60) {phDiff=72-phDiff;}
        //
        etMargin = max((etResult.mag/100), (unsigned) 1) + 2;
        if (etResult.mag==0) { phMargin = 72; } else { phMargin = (30/etResult.mag) + 1; }
        if ((etDiff > etMargin) || (phDiff > phMargin)) {cout << "Algo etMiss diff "
                                                    << etDiff << " phi diff " << phDiff << endl; testPass = false;}
      }
      // Check the other output values
      if (!myGlobalEnergy->getEtSum().overFlow() && !etTotalOvrFlow &&
          (myGlobalEnergy->getEtSum().value() != etTotal)) {cout << "Algo etSum" << endl; testPass = false;}
      if (!myGlobalEnergy->getEtHad().overFlow() && !htTotalOvrFlow &&
          (myGlobalEnergy->getEtHad().value() != htTotal)) {cout << "Algo etHad" << endl; testPass = false;}
//       for (int j=0 ; j<12 ; j++) {
//         if (myGlobalEnergy->getJetCount(j).value() != JcResult[j]) {cout << "Algo jCount" << endl; 
//       cout << "Expected " << JcResult[j] << " found " << myGlobalEnergy->getJetCount(j) << endl; 
//       cout << "PlusWheel " << myGlobalEnergy->getInputJcValPlusWheel(j) << endl; 
//       cout << "PlusWheel " << myGlobalEnergy->getInputJcVlMinusWheel(j) << endl; testPass = false;}
//       }
    
      //
      //--------------------------------------------------------------------------------------
      //
    
      if(testPass == false)
	{
	  // Print failed events for debug purposes
	  cout << *myGlobalEnergy << endl;
	  cout << "Failed at test number " << t <<endl;
	  return "Test class has failed algorithm processing!";
	}
    }
  // Test the print routine at the end
  cout << *myGlobalEnergy << endl;
  return "Test class has passed!";
}

// Generates 2-d missing Et vector
void generateMissingEtTestData(int &Ex, int &Ey, etmiss_vec &Et)
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
  unsigned dummySum;
  float p,r,s;
  unsigned Emag, Ephi;

  const float nbins = 18.;

  // Generate a pair of uniform pseudo-random integers
  generateTestData(components, (int) 2, rmax, dummySum);

  // Exclude the value zero for the first random integer
  // (Alternatively, return an overflow bit)
  while (components[0]==0) {generateTestData(components, (int) 2, rmax, dummySum);}

  // Convert to the 2-d Gaussian
  r = float(rmax);
  s = r/float(components[0]);
  p = float(components[1])/r;
  // Force phi value into the centre of a bin
  Emag = int(sigma*sqrt(2.*log(s)));
  Ephi = int(nbins*p);
  //
  Et.mag = Emag;
  Et.phi = (4*Ephi)+2;

  Ex = etComponent(Emag, ((2*Ephi)+10)%36);
  Ey = etComponent(Emag, ((2*Ephi)+1));
}

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

/// Check the initial Ht calculation
bool checkHt(L1GctJetLeafCard* jlc, unsigned &leafHt) {

  unsigned result=jlc->getOutputHt().value();
  unsigned sumHt=0;

  for (int i=0; i<L1GctJetFinder::MAX_JETS_OUT; i++) {
    if (!jlc->getOutputJetsA().at(i).isNullJet()) { sumHt += jetHtSum(jlc->getJetFinderA(), i); }
    if (!jlc->getOutputJetsB().at(i).isNullJet()) { sumHt += jetHtSum(jlc->getJetFinderB(), i); }
    if (!jlc->getOutputJetsC().at(i).isNullJet()) { sumHt += jetHtSum(jlc->getJetFinderC(), i); }
  }

  leafHt = result;
  return (sumHt==result);
}

/// Work out the Ht for a jet.
unsigned jetHtSum(L1GctJetFinder* jf, int jn) {

  // We need to take the eta and phi, and go back to the
  // raw data, since the raw energy sum is not stored
  vector<L1GctRegion>inputRegions = jf->getInputRegions();
  const unsigned COL_OFFSET = ((L1GctMap::N_RGN_ETA)/2)+1;

  // Check the input array size
  if (inputRegions.size()!=(COL_OFFSET*4)) {
    cout << "Invalid size for jet finder input " << inputRegions.size() << " expecting " << (COL_OFFSET*4) << endl;
    return 0;
  }

  // Find the eta and phi for this jet, and check
  unsigned eta = static_cast<unsigned>(jf->getJets().at(jn).eta());
  unsigned phi = static_cast<unsigned>(jf->getJets().at(jn).phi());

  if (phi>1 || eta>=(COL_OFFSET-1)) {
    cout << "Invalid eta, phi for jet: " << eta << ", " << phi << endl;
    return 0;
  }

  // Sum the et values for the nine regions centred on this eta and phi
  unsigned rawSum = 0;

  for (unsigned col=phi; col<=phi+2; col++) {
    rawSum += inputRegions.at(col*COL_OFFSET+eta).et();
    rawSum += inputRegions.at(col*COL_OFFSET+eta+1).et();
    if ((eta+2)<COL_OFFSET) {
      rawSum += inputRegions.at(col*COL_OFFSET+eta+2).et();
    }
  }
  // Convert to Ht and return
  return static_cast<unsigned>(jf->getJetEtCalLut()->convertToTenBitRank(static_cast<uint16_t>(rawSum),
									 static_cast<uint16_t>(eta)));
}


/// Generates test data consisting of energies to be added together with their sum
void generateTestData(vector<unsigned> &energies, int size, unsigned max, unsigned &sum)
{
  int i;
  int r,e;
  float p,q,s,t;
  int tempsum;

  tempsum = 0;
  p = float(max);
  q = float(RAND_MAX);
  for (i=0; i<size; i++) {
    r = rand();
    s = float(r);
    t = s*p/q;
    e = int(t);

    energies[i] = e;
    tempsum = tempsum+e;
  }
  sum = tempsum;
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
