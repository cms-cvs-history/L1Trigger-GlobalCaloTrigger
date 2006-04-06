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

// //Custom headers needed for this test
// #include "L1Trigger/GlobalCaloTrigger/interface/L1GctRegion.h"
// #include "L1Trigger/GlobalCaloTrigger/interface/L1GctJet.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <exception> //for exception handling
#include <stdexcept> //for std::runtime_error()
using namespace std;

// //Typedefs for the vector templates used
// typedef vector<L1GctRegion> RegionsVector;
// typedef vector<L1GctJet> JetsVector;

//  FUNCTION PROTOTYPES
/// Runs the test, and returns a string with the test result message in.
string classTest();
//
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
  L1GctGlobalEnergyAlgos myGlobalEnergy;
  bool testPass = true;       //Test passing flag.
    
  const int noOfEtValues=2;
  const int noOfHtValues=3;
  const int maxValues=3;
  vector<unsigned> energyValues(maxValues);
  unsigned energySum;

  unsigned EtSumResult;
  unsigned EtHadResult;
  vector<unsigned> JcResult(12);

  const int nbitsEt=12;
  const unsigned energyMax=((1<<nbitsEt) - 1);

  const int maxTests=1000;

  for (int t=0; t<maxTests; t++)
    {
      // Initialise the myGlobalEnergy object
      myGlobalEnergy.reset();

      // For each of the different types of energy sum,
      // generate our test input data and known results.
      // Load the input data into the algorithm and store
      // a local copy of the output.
      //
      //--------------------------------------------------------------------------------------
      //
      // Total transverse energy (Et):
      generateTestData(energyValues, noOfEtValues, energyMax, energySum);

      // Fill the L1GctGlobalEnergyAlgos
      myGlobalEnergy.setInputWheelEt(0, energyValues[0], false);
      myGlobalEnergy.setInputWheelEt(1, energyValues[1], false);

      // Test the GetInput... methods
      if (myGlobalEnergy.getInputEtValPlusWheel() != energyValues[0]) { testPass = false; }
      if (myGlobalEnergy.getInputEtVlMinusWheel() != energyValues[1]) { testPass = false; }
            
      // Local storage of the sum
      if (energySum <= energyMax)
	{ EtSumResult = energySum; }
      else
	{ EtSumResult = (energySum & energyMax) | (1<<nbitsEt); }

      //
      //--------------------------------------------------------------------------------------
      //
      // Jet energy sums (Ht):
      generateTestData(energyValues, noOfHtValues, energyMax, energySum);

      // Fill the L1GctGlobalEnergyAlgos
      myGlobalEnergy.setInputWheelHt(0, energyValues[0], false);
      myGlobalEnergy.setInputWheelHt(1, energyValues[1], false);
      myGlobalEnergy.setInputBoundaryHt(energyValues[2], false);

      // Test the GetInput... methods
      if (myGlobalEnergy.getInputHtValPlusWheel() != energyValues[0]) { testPass = false; }
      if (myGlobalEnergy.getInputHtVlMinusWheel() != energyValues[1]) { testPass = false; }
      if (myGlobalEnergy.getInputHtBoundaryJets() != energyValues[2]) { testPass = false; }
            
      // Local storage of the sum
      if (energySum <= energyMax)
	{ EtHadResult = energySum; }
      else
	{ EtHadResult = (energySum & energyMax) | (1<<nbitsEt); }

      //
      //--------------------------------------------------------------------------------------
      //
      // Jet counts
      for (int j=0; j<12; j++) {
	unsigned total;
        vector <unsigned> values(3);
	const unsigned maxvalue=20;

        generateTestData(values, (int) 3, maxvalue, total);

        myGlobalEnergy.setInputWheelJc(0, j, values[0]);
	myGlobalEnergy.setInputWheelJc(1, j, values[1]);
        myGlobalEnergy.setInputBoundaryJc(j, values[2]);

        if (values[0]>=15) {
	  values[0] = 15;
	  total = 31;
	}
        if (values[1]>=15) {
	  values[1] = 15;
	  total = 31;
	}
        if (values[2]>=7) {
	  values[2] = 7;
	  total = 31;
	}
        if (total>=32) total = 31;
 
        if (myGlobalEnergy.getInputJcValPlusWheel(j) != values[0]) {testPass = false;}
        if (myGlobalEnergy.getInputJcVlMinusWheel(j) != values[1]) {testPass = false;}
        if (myGlobalEnergy.getInputJcBoundaryJets(j) != values[2]) {testPass = false;}

        JcResult[j] = total;
      }

      if(testPass == false)
	{
	  return "Test class has failed initial data input/output comparison!";
	}
      //
      //--------------------------------------------------------------------------------------
      //
   
      myGlobalEnergy.process();  //Run algorithm

      //
      //--------------------------------------------------------------------------------------
      //
      // Check the output values
      if (myGlobalEnergy.getEtSum() != EtSumResult) {testPass = false;}
      if (myGlobalEnergy.getEtHad() != EtHadResult) {testPass = false;}
      for (int j=0 ; j<12 ; j++) {
        if (myGlobalEnergy.getJetCounts(j) != JcResult[j]) {testPass = false;}
      }
    
      //
      //--------------------------------------------------------------------------------------
      //
    
      if(testPass == false)
	{
	  return "Test class has failed algorithm processing!";
	}
    }
    return "Test class has passed!";        
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