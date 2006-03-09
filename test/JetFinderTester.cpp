/*! \file JetFinderTester.cpp
 * \brief Procedural skeleton unit-test code for the L1GctJetFinder class.
 *
 *  This is skeleton code that reads in data from a text file to feed into
 *  the setInputRegions() method, runs the process() method, and then
 *  checks the data from the outputting methods getInputRegions() and
 *  getJets() against known results also stored in the text file.
 *
 * \author Robert Frazier
 * \date March 2006
 */


//NOTE all these includes need to be sorted with a proper build include path,
//rather than a relative path, in order to comply with CMS style.

#include "../interface/L1GctJetFinder.h"  //The class to be tested

//Custom headers needed for this test
#include "../interface/L1GctRegion.h"
#include "../interface/L1GctJet.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
using namespace std;

//Typedefs for the vector templates used
typedef vector<L1GctRegion> RegionsVector;
typedef vector<L1GctJet> JetsVector;

//  FUNCTION PROTOTYPES
/// Loads test input regions and also the known results from a text file.
bool loadTestData(RegionsVector& regions, JetsVector& jets, ifstream &fin, const string &fileName);
/// Function to safely open files of any name, using a referenced return ifstream
bool safeOpenFile(ifstream &fin, const string &name);

/// Start of unit test code
int main(int argc, char **argv)
{
    cout << "\nL1GctJetFinder class unit tester." << endl;

    L1GctJetFinder myJetFinder; //TEST OBJECT;    
    bool testPass = true;       //Test passing flag.
    
    // Number of calorimter regions to be fed to the jet finder.
    const int maxRegions = 64;  //64 based on old GCT design
    
    // Name of the file containing the test input data.
    const string testDataFile = "JetFinderTesterData.txt";  
    
    // File input stream
    ifstream fin;
    
    // Vectors for reading in test data from the text file.
    RegionsVector inputRegions(maxRegions);
    JetsVector correctJets;
    // Vectors for receiving the output from the object under test.
    RegionsVector outputRegions(maxRegions);
    JetsVector outputJets;
    
    // Load our test input data and known results
    if(!loadTestData(inputRegions, correctJets, fin, testDataFile)) { return 0; }
    
    //Fill the L1GctJetFinder with regions.
    for(int i = 0; i < maxRegions; ++i)
    {
        myJetFinder.setInputRegion(i, inputRegions[i]);
    }

    // Test the getInputRegion method
    outputRegions = myJetFinder.getInputRegions();
    for(int i = 0; i < maxRegions; ++i)
    {
        if(outputRegions[i].getEt() != inputRegions[i].getEt()) { testPass = false; break; }
        if(outputRegions[i].getMip() != inputRegions[i].getMip()) { testPass = false; break; }
        if(outputRegions[i].getQuiet() != inputRegions[i].getQuiet()) {testPass = false; break; }
    }
    
    if(testPass == false)
    {
        cout << "\nTest class has failed initial data input/output!" << endl;
        return 0;
    }
    
    myJetFinder.process();  //Run algorithm
    
    //Get and then test the output jets against known results
    //NEEDS FILLING IN SIMILARLY TO ABOVE
    outputJets = myJetFinder.getJets();
    // **Fill in***
    
    
    if(testPass == false)
    {
        cout << "\nTest class has failed algorithm processing!" << endl;
        return 0;
    }    
        
    cout << "\nTest class has passed!" << endl;   
    return 0;   
}


// Loads test input regions from a text file.
bool loadTestData(RegionsVector& regions, JetsVector& jets, ifstream &fin, const string &fileName)
{
    if(!safeOpenFile(fin, fileName)) {return false;}  //open the file or fail gracefully
    
    unsigned long int tempEt = 0;
    unsigned short int tempMip = false;
    unsigned short int tempQuiet = false;
    
    //loads the input data
    for(unsigned int i = 0; i < regions.size(); ++i)
    {
        //read in the data from the file
        fin >> tempEt;
        fin >> tempMip;
        fin >> tempQuiet;
        
        regions[i].setEt(tempEt);
        if(tempMip == 0) { regions[i].setMip(false); } else { regions[i].setMip(true); }
        if(tempQuiet == 0) { regions[i].setQuiet(false); } else { regions[i].setQuiet(true); }
    }
    
    //do similar to load the 'known' output jets (that we currently don't know...)
    
        
    return true;
}
    
    
    
// Function to safely open files of any name, using a referenced return ifstream
bool safeOpenFile(ifstream &fin, const string &name)
{
  //Opens the file
  fin.open(name.c_str(), ios::in);

  //Error message, and return false if it goes pair shaped
  if(!fin.good())
  {
    cout << "\nError! Couldn't open the file "<< name << endl;
    return false;
  }
  return true;
}
