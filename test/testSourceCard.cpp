/*! \file testSourceCard.cpp
 * \brief Partial test code for testing L1GctSourceCard
 *
 *  Reads in an RCT output file using the three variants of an
 *  L1GctSourceCard, then outputs the data back out to file
 *  in the same format.  This isn't full unit-test code.
 *
 * \author Robert Frazier
 * \date April 2006
 */

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <fstream>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <bitset>
using namespace std;

typedef vector<L1GctEmCand> ElecVector;
typedef vector<L1GctRegion> Region;
typedef unsigned long int ULong;

//Global vectors, bitsets, etc, to hold output data
ULong currentBX;
vector<L1GctEmCand> isoElectronsOut;  
vector<L1GctEmCand> nonIsoElectronsOut;
unsigned mipBitsOut;
unsigned quietBitsOut;
vector<L1GctRegion> regionsOutT2;
vector<L1GctRegion> regionsOutT3;

//Global file output
ofstream fout;

// Name of the files for test input data and results output.
const string testDataFile = "testSourceCardInput.txt";  
const string resultsFile = "testSourceCardOutput.txt";


//  FUNCTION PROTOTYPES
/// Outputs to file in RCT output format the data that the 3 source card variants have read in
void outputDataToFile(ofstream& fout);
/// Function to safely open output files of any name, using a referenced return ofstream
void safeOpenOutputFile(ofstream &fout, const string name);
/// Converts the EmCand info back into the original RCT output format and returns as a ULong
ULong compressElectronData(L1GctEmCand& emCand);
/// Converts central Region info back into the original RCT output format and returns as a ULong
ULong compressCentralRegionData(L1GctRegion& region);


int main(int argc, char **argv)
{
  cout << "\n********************************************" << endl;
  cout << "L1GctSourceCard class (partial) unit tester." << endl;
  cout << "********************************************" << endl;    
  
  try
  {
    safeOpenOutputFile(fout, resultsFile);
    
    L1GctSourceCard * mySourceCardT1 = new L1GctSourceCard(0, L1GctSourceCard::cardType1); //test objects
    L1GctSourceCard * mySourceCardT2 = new L1GctSourceCard(0, L1GctSourceCard::cardType2); 
    L1GctSourceCard * mySourceCardT3 = new L1GctSourceCard(0, L1GctSourceCard::cardType3);
      
    //cout << mySourceCardT1 << endl;
    //cout << mySourceCardT2 << endl;
    //cout << mySourceCardT3 << endl;        

    //Read the same file simultaneously
    mySourceCardT1->openInputFile(testDataFile);
    mySourceCardT2->openInputFile(testDataFile);  
    mySourceCardT3->openInputFile(testDataFile);        
    
    while(mySourceCardT1->dataValid())
    {        
      mySourceCardT1->readBX(); //Get the data for each bunch crossing for each card type
      mySourceCardT2->readBX();
      mySourceCardT3->readBX();
      
      currentBX = mySourceCardT1->getBxNum();
      
      //Read out the data
      isoElectronsOut = mySourceCardT1->getIsoElectrons();
      nonIsoElectronsOut = mySourceCardT1->getNonIsoElectrons();
      mipBitsOut = mySourceCardT1->getMipBits();
      quietBitsOut = mySourceCardT1->getQuietBits();
      regionsOutT2 = mySourceCardT2->getRegions();
      regionsOutT3 = mySourceCardT3->getRegions();           
      
      //output combined data to file
      outputDataToFile(fout);
    }
    
    fout.close();
    
    delete mySourceCardT1;
    delete mySourceCardT2;
    delete mySourceCardT3;
    
    cout << "See output file " << resultsFile << endl;
  }
  catch (cms::Exception& e)
  {
    cerr << e.what() << endl;
  }
  catch(...)
  {
    cerr << "\nError! An unknown exception has occurred!" << endl;
  }
  
  return 0;   
}

// Outputs to file in RCT output format the data that the 3 source card variants have read in
void outputDataToFile(ofstream& fout)
{
  unsigned int i=0; //counter
  dec(fout);
  fout << "Crossing  " << currentBX << endl;
  hex(fout);
  for(i=0;i<4;++i)
  {
    fout << compressElectronData(isoElectronsOut[i]) << " ";
  }
  for(i=0;i<4;++i)
  {
    fout << compressElectronData(nonIsoElectronsOut[i]) << " ";
  }
  fout << endl;
  for(i=0;i<14;++i)
  {
    fout << " " << ((mipBitsOut & (1<<i)) >> i);
  }
  fout << endl;
  for(i=0;i<14;++i)
  {
    fout << " " << ((quietBitsOut& (1<<i)) >> i);
  }
  fout << endl;   
  for(i=0;i<regionsOutT3.size();++i)
  {
    fout << " " << compressCentralRegionData(regionsOutT3[i]);
  }
  for(i=0;i<regionsOutT2.size();++i)
  {
    if(i==4) { fout << endl; }
    
    if(i<4) { fout << " " << compressCentralRegionData(regionsOutT2[i]); }
    else    { fout << " " << regionsOutT2[i].et(); }
  }
  fout << endl;
}


// Function to safely open output files of any name, using a referenced return ofstream
void safeOpenOutputFile(ofstream &fout, const string name)
{
  //Opens the file
  fout.open(name.c_str(), ios::trunc);
  
  //Throw an exception if something is wrong
  if(!fout.good())
  {
    throw cms::Exception("FileWriteError")
    << "Couldn't open the file " << name << " for writing!\n";
  }
  return;
}

// Converts the EmCand info back into the original RCT output format and returns as a ULong
ULong compressElectronData(L1GctEmCand& emCand)
{
  ULong rank = emCand.rank();
  ULong phi  = emCand.phi();
  ULong eta  = emCand.eta();
  
  phi <<= 6;
  eta <<= 7;
  
  return rank+phi+eta;
}

// Converts central Region info back into the original RCT output format and returns as a ULong
ULong compressCentralRegionData(L1GctRegion& region)
{
  ULong output=0;
  if(region.overFlow()) {output += (1<<10);}
  if(region.tauVeto()) {output += (1<<11);}
  output += region.et();
  
  return output;    
}
