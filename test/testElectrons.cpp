/*! \file testElectrons.cpp
 * \test file for testing the full GCT chain for electrons. 
 *
 * 
 * \author Maria Hansen
 * \date April 2006
 */


#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctElectronSorter.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEmLeafCard.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "FWCore/Utilities/interface/Exception.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <sstream>  //for int->char conversion
#include <exception>
#include <stdexcept>

using namespace std;

//  Function for reading in the dummy data                                
void LoadFileData(const string &inputFile);
//Function to easily output a gctEmCand vector to screen                              
void print(vector<L1GctEmCand> cands);

vector<L1CaloEmCand> data;
vector<L1GctEmCand> gctData;
vector<L1CaloEmCand> fileIso;
vector<L1CaloEmCand> fileNoniso;
ifstream file;


int main()
{
  cout <<"\n=============================="<<endl;
  cout <<" Test full chain of electrons "<<endl;
  cout <<"=============================="<<endl;

  vector<L1GctEmCand> isoElectrons;
  vector<L1GctEmCand> nonIsoElec;
  vector<L1GctEmLeafCard*> emLeafs;
  bool isos = false;
  bool nonisos = false;
  bool leafs = false;

  try { 
    L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger();
    L1GctElectronSorter* isoSort = new L1GctElectronSorter(72,0);
    L1GctElectronSorter* nonIsoSort = new L1GctElectronSorter(72,1);

    isoElectrons = gct->getIsoElectrons();
    nonIsoElec = gct->getNonIsoElectrons();
    emLeafs = gct->getEmLeafCards();

    // Check that internal vectors in GCT is correctly sat up for electrns
    if(isoElectrons.size()!= 4){
      throw cms::Exception("ErrorSizeOfIso")
	<< "The GCT class returns the wrong number of iso electrons"<<endl;
      isos = true;
    }
    if(nonIsoElec.size()!= 4){
      throw cms::Exception("ErrorSizeOfNonIso")
	<< "The GCT class returns the wrong number of non-iso electrons"<<endl;
      nonisos = true;
    }
    if(emLeafs.size()!= 2){
      throw cms::Exception("ErrorSizeOfEmLeaf")
	<< "The GCT class holds the wrong number of leaf cards"<<endl;
      leafs = true;
    }
  
    //Open source card files and run through the gct chain, returning the four electrons highest in rank
    // testElectronsRct_x are random test files with rct data
    //testElectrons_x are testfiles containing 1,2,3 electrons, equal rank,phi and energies etc.
    // testEmCandsDummy_ are empty test files 

    //    std::string fileName = "data/testElectronsRct_";
    std::string fileName = "data/testEmDummy_";

    gct->openSourceCardFiles(fileName);
    gct->process();
    vector<L1GctEmCand> newIso = gct->getIsoElectrons();
    vector<L1GctEmCand> newnonIso = gct->getNonIsoElectrons();
    
    cout<<"=========== From the GCT chain ==============="<<endl;
    cout<<"Iso electrons are: "<<endl;
    print(newIso);
    cout<<"Non-iso electrons are: "<<endl;
    print(newnonIso);

    //Open the same files using the function LoadFile and sort them and see if same output it returned
    for(int i=0;i<18;i++){
      stringstream ss;
      string fileNo;
      ss << i;
      ss >> fileNo;
      LoadFileData(fileName+fileNo);
    }
    for(int unsigned n=0;n!=18;n++){
      fileIso.push_back(data[n*8]);
      fileIso.push_back(data[n*8 + 1]);
      fileIso.push_back(data[n*8 + 2]);
      fileIso.push_back(data[n*8 + 3]);
      fileNoniso.push_back(data[n*8 + 4]);
      fileNoniso.push_back(data[n*8 + 5]);
      fileNoniso.push_back(data[n*8 + 6]);
      fileNoniso.push_back(data[n*8 + 7]);
    }
 
    for(unsigned int i=0;i!=72;i++){
      isoSort->setInputEmCand(i,fileIso[i]);
      nonIsoSort->setInputEmCand(i,fileNoniso[i]);
    }
    isoSort->process();
    nonIsoSort->process();
    vector<L1GctEmCand> iso = isoSort->getOutputCands();
    vector<L1GctEmCand> noniso = nonIsoSort->getOutputCands();
    cout<<"=================From files externally sorted============="<<endl;
    cout<<"Iso electrons are:"<<endl;
    print(iso);
    cout<<"Non-iso electrons are:"<<endl;
    print(noniso);
   
    if(!isos && !nonisos && !leafs){
      cout<<"======================================="<<endl;
      cout<<" Test full chain of electrons ran"<<endl;
      cout<<" without errors"<<endl;
      cout<<"========================================"<<endl;
    }

    delete isoSort;
    delete nonIsoSort;
    delete gct;

  }
  catch(std::exception e){
    cerr << e.what() << endl;
  }
  return 0;
}

// Function definition of function that reads in dummy data and load it into inputCands vector        
void LoadFileData(const string &inputFile)
{
  //Opens the file                                                                                     
  file.open(inputFile.c_str(), ios::in);

  if(!file){
    throw cms::Exception("ErrorOpenFile")
      << "Cannot open input data file" << endl;
  }

  unsigned candRank = 0, candRegion = 0, candCard = 0, candCrate = 0;
  short dummy;
  string bxNo = "poels";
  bool candIso=0;

  //Reads in first crossing, then crossing no                                           
  //then 8 electrons, 4 first is iso, next non iso.                                                    
  //The 3x14 jet stuff and 8 mip and quiet bits are skipped over.                                     
   
    file>> bxNo;
    file >>std::dec>> dummy;
    for(int i=0; i<58; i++){ //Loops over one bunch-crossing                                          
      if(i<8){
        if(i>3){
          file >>std::hex>> dummy;
          candRank = dummy & 0x3f;
          candRegion = (dummy>>6) & 0x1;
          candCard = (dummy>>7) & 0x7;
          candIso = 1;
          L1CaloEmCand electrons(candRank, candRegion, candCard, candCrate, candIso);
          data.push_back(electrons);
	}else{
          if(i<4){
            file >>std::hex>> dummy;
            candRank = dummy & 0x3f;
            candRegion = (dummy>>6) & 0x1;
            candCard = (dummy>>7) & 0x7;
            candIso = 0;
            L1CaloEmCand electrons(candRank, candRegion, candCard, candCrate, candIso);
            data.push_back(electrons);
          }else{
	    file>>dummy;
          }
        }
      }else{
        file >>std::hex>> dummy;
      }
    }
  file.close();
  return;
}

void print(vector<L1GctEmCand> cands){
  for(unsigned int i=0; i!=cands.size(); i++){
    cout<<"          Rank: "<<cands[i].rank()<<"  Eta: "<<cands[i].etaIndex()<<"  Phi: "<<cands[i].phiIndex()<<endl;
  }
  return;
}


