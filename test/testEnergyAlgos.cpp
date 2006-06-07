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
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h"

// //Custom headers needed for this test
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctRegion.h"
// #include "L1Trigger/GlobalCaloTrigger/interface/L1GctJet.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEtTypes.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <exception> //for exception handling
#include <stdexcept> //for std::runtime_error()
using namespace std;

// //Typedefs for the vector templates used
        struct etmiss_vec { unsigned mag; unsigned phi;};
typedef vector<L1GctRegion> RegionsVector;
// typedef vector<L1GctJet> JetsVector;

//  FUNCTION PROTOTYPES
/// Runs the test, and returns a string with the test result message in.
string classTest();
//
/// Generates test data for missing Et both as (x,y) components and
/// 2-vector (magnitude, direction)
void generateMissingEtTestData(int &Ex, int &Ey, etmiss_vec &Et);
/// Generates test data consisting of energies to be added together with their sum
void generateTestData(vector<unsigned> &energies, int size, unsigned max, unsigned &sum);

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
  vector<L1GctSourceCard*> mySourceCards;

  bool testPass = true;       //Test passing flag.
    
  RegionsVector inputRegions;
  const int noOfInputRegions=10;

//   const int noOfEtValues=2;
//   const int noOfHtValues=2;
  const int maxValues=3;
  vector<unsigned> energyValues(maxValues);
//   unsigned energySum;

//   unsigned EtSumResult;
//   unsigned EtHadResult;
  vector<unsigned> JcResult(12);

//   const int nbitsEt=12;
//   const unsigned energyMax=((1<<nbitsEt) - 1);

  int exPlusVal, exMinusVl;
  int eyPlusVal, eyMinusVl;
  unsigned etPlusVal, etMinusVl;
  int exValue, eyValue;
  etmiss_vec etVector;
  bool etOvflo;
  unsigned etDiff, phDiff;
  unsigned etMargin, phMargin;

  const int maxTests=100000;
  const int initialTests=10;

  // Initialise the gct
  L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(false);
  myGlobalEnergy = gct->getEnergyFinalStage();
  mySourceCards  = gct->getSourceCards();

  // Initialise my input regions
  for (int i=0; i<noOfInputRegions; i++) {
    inputRegions.push_back(L1GctRegion(0, 0, 0, false, false, false, false));
  }

  for (int t=0; t<maxTests; t++)
    {
      if (t<initialTests) { cout << "Test " << t << endl; }
      // Initialise the gct
      gct->reset();

      exPlusVal = 0;
      eyPlusVal = 0;
      etPlusVal = 0;
      exMinusVl = 0;
      eyMinusVl = 0;
      etMinusVl = 0;

      for (int i=0; i<(t<initialTests ? 1 : 0); i++) {

        // Initial test: Set non-zero energy in a single input region and check
        // that it is the same at the output (within rounding errors)
        // Later tests (to be added): repeat the check for multiple regions
        //
        //--------------------------------------------------------------------------------------
        //
        // Missing Et (Ex, Ey):
        generateMissingEtTestData(exValue, eyValue, etVector);
	cout << "Region energy " << etVector.mag << " angle bin " << etVector.phi << endl;
	cout << "Components ex " << exValue << " ey " << eyValue << endl;
        //
        // Fill the Source Card input
        // Set a single region input
        unsigned phiRegion = etVector.phi/4;
        L1GctRegion temp(0, phiRegion%2, etVector.mag, false, false, false, false);
        if (temp.phi()==0) {
	  inputRegions[0] = temp;
	  inputRegions[6] = 0;
        } else {
	  inputRegions[0] = 0;
	  inputRegions[6] = temp;
        }
        mySourceCards[(3*(phiRegion/2))+2]->setRegions(inputRegions);

        // Here we fill the expected values
        exMinusVl += exValue;
        eyMinusVl += eyValue;
        etMinusVl += etVector.mag;
      }

      //
      // Run the processing
      //--------------------------------------------------------------------------------------
      //
   
      gct->process();  //Run algorithm

      //
      // Check the input to the final GlobalEnergyAlgos is as expected
      //--------------------------------------------------------------------------------------
      //
      if (myGlobalEnergy->getInputExVlMinusWheel().value()!=exMinusVl) { testPass = false; }
      if (myGlobalEnergy->getInputExValPlusWheel().value()!=exPlusVal) { testPass = false; }
      if (myGlobalEnergy->getInputEyVlMinusWheel().value()!=eyMinusVl) { testPass = false; }
      if (myGlobalEnergy->getInputEyValPlusWheel().value()!=eyPlusVal) { testPass = false; }
      if (myGlobalEnergy->getInputEtVlMinusWheel().value()!=etMinusVl) { testPass = false; }
      if (myGlobalEnergy->getInputEtValPlusWheel().value()!=etPlusVal) { testPass = false; }
   
      if(testPass == false)
	{
	  cout << *myGlobalEnergy << endl;
	  return "Test class has failed initial data input/output comparison!";
	}
      //
      //--------------------------------------------------------------------------------------
      //
   
      cout << "Ex, Ey input values" << endl;
      cout << "Minus wheel Ex 0 " << gct->getWheelEnergyFpgas()[0]->getInputEx(0) << endl;
      cout << "Minus wheel Ex 1 " << gct->getWheelEnergyFpgas()[0]->getInputEx(1) << endl;
      cout << "Minus wheel Ex 2 " << gct->getWheelEnergyFpgas()[0]->getInputEx(2) << endl;
      cout << " Plus wheel Ex 0 " << gct->getWheelEnergyFpgas()[1]->getInputEx(0) << endl;
      cout << " Plus wheel Ex 1 " << gct->getWheelEnergyFpgas()[1]->getInputEx(1) << endl;
      cout << " Plus wheel Ex 2 " << gct->getWheelEnergyFpgas()[1]->getInputEx(2) << endl;
      cout << "Minus wheel Ey 0 " << gct->getWheelEnergyFpgas()[0]->getInputEy(0) << endl;
      cout << "Minus wheel Ey 1 " << gct->getWheelEnergyFpgas()[0]->getInputEy(1) << endl;
      cout << "Minus wheel Ey 2 " << gct->getWheelEnergyFpgas()[0]->getInputEy(2) << endl;
      cout << " Plus wheel Ey 0 " << gct->getWheelEnergyFpgas()[1]->getInputEy(0) << endl;
      cout << " Plus wheel Ey 1 " << gct->getWheelEnergyFpgas()[1]->getInputEy(1) << endl;
      cout << " Plus wheel Ey 2 " << gct->getWheelEnergyFpgas()[1]->getInputEy(2) << endl;
      cout << "Minus wheel Et 0 " << gct->getWheelEnergyFpgas()[0]->getInputEt(0) << endl;
      cout << "Minus wheel Et 1 " << gct->getWheelEnergyFpgas()[0]->getInputEt(1) << endl;
      cout << "Minus wheel Et 2 " << gct->getWheelEnergyFpgas()[0]->getInputEt(2) << endl;
      cout << " Plus wheel Et 0 " << gct->getWheelEnergyFpgas()[1]->getInputEt(0) << endl;
      cout << " Plus wheel Et 1 " << gct->getWheelEnergyFpgas()[1]->getInputEt(1) << endl;
      cout << " Plus wheel Et 2 " << gct->getWheelEnergyFpgas()[1]->getInputEt(2) << endl;
      cout << endl;
      cout << "Minus wheel Ex total " << myGlobalEnergy->getInputExVlMinusWheel() << endl;
      cout << " Plus wheel Ex total " << myGlobalEnergy->getInputExValPlusWheel() << endl;
      cout << "Minus wheel Ey total " << myGlobalEnergy->getInputEyVlMinusWheel() << endl;
      cout << " Plus wheel Ey total " << myGlobalEnergy->getInputEyValPlusWheel() << endl;
      cout << "Minus wheel Et total " << myGlobalEnergy->getInputEtVlMinusWheel() << endl;
      cout << " Plus wheel Et total " << myGlobalEnergy->getInputEtValPlusWheel() << endl;
      cout << endl;
      cout << "Output Et Miss       " << myGlobalEnergy->getEtMiss() << endl;
      cout << "Output Et Miss angle " << myGlobalEnergy->getEtMissPhi() << endl;
      cout << "Output Et Total      " << myGlobalEnergy->getEtSum() << endl;
      cout << endl;
      cout << "Values loaded magnitude " << etVector.mag << " angle bin " << etVector.phi << endl;
      cout << endl;

      //
      //--------------------------------------------------------------------------------------
      //
      // Check the missing Et calculation. Allow some margin for the
      // integer calculation of missing Et.
      if (etOvflo) {
	if (!myGlobalEnergy->getEtMiss().overFlow()) {testPass = false;}
      } else {
        etDiff = (unsigned) abs((long int) etVector.mag - (long int) myGlobalEnergy->getEtMiss().value());
        phDiff = (unsigned) abs((long int) etVector.phi - (long int) myGlobalEnergy->getEtMissPhi().value());
        if (phDiff>60) {phDiff=72-phDiff;}
        //
        etMargin = max((etVector.mag/100), (unsigned) 1) + 5;
        phMargin = (30/etVector.mag) + 3;
        if ((etDiff > etMargin) || (phDiff > phMargin)) {cout << "Algo etMiss diff "
                                                    << etDiff << " phi diff " << phDiff << endl; testPass = false;}
      }
//       // Check the other output values
//       if (myGlobalEnergy->getEtSum().value() != EtSumResult) {cout << "Algo etSum" << endl; testPass = false;}
//       if (myGlobalEnergy->getEtHad().value() != EtHadResult) {cout << "Algo etHad" << endl; 
//       cout << "Expected " << EtHadResult << " found " << myGlobalEnergy->getEtHad() << endl; testPass = false;}
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
// 	  cout << *myGlobalEnergy << endl;
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
  const float sigma=400.;

  // rmax controls the magnitude range
  // Chosen as a power of two conveniently close to
  // exp(5*5/2) to give a 5*sigma range.
  const unsigned rmax=262144;

  vector<unsigned> components(2);
  unsigned dummySum;
  float p,r,s;
  float Emag, Ephi;

  // Generate a pair of uniform pseudo-random integers
  generateTestData(components, (int) 2, rmax, dummySum);

  // Exclude the value zero for the first random integer
  // (Alternatively, return an overflow bit)
  while (components[0]==0) {generateTestData(components, (int) 2, rmax, dummySum);}

  // Convert to the 2-d Gaussian
  r = float(rmax);
  s = r/float(components[0]);
  p = float(components[1])/r;
  Emag = sigma*sqrt(2.*log(s));
  Ephi = 6.2831854*p;
  //
  Et.mag = int (Emag);
  Et.phi = int (72.*p);
  Ex = int (Emag*cos(Ephi));
  Ey = int (Emag*sin(Ephi));
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
