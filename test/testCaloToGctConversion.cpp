
/*! \file tesCaloToGctConversion.cpp
 * \test file for testing the conversion between calo and gct candidates
 * 
 *  This test program reads in dummy data, followed by testing
 *  its methods are working correctly. Any discrepancies found
 *  between read in data and output will be written to the screen
 *  
 *
 * \author Maria Hansen
 * \date July 2006
 */

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctElectronSorter.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <exception>
#include <stdexcept>

using namespace std;

// typedefs and others
typedef vector<L1CaloEmCand> EmInputCandVec;
typedef vector<L1GctEmCand> EmOutputCandVec;
ifstream file;
ofstream ofile;
vector<L1CaloEmCand> data;
vector<L1GctEmCand> gctData;
// Function that prints out the eta and phi value for a given (crateNo) Rct crate
void eta_phiTable(int crateNo);
//  Function for reading in the dummy data into the caloEmcand vector data
void LoadFileData(const string &inputFile, bool iso);
// Function to print out gct candidate data to the screen
void print(EmOutputCandVec toPrint);
// Copy of function in electron sorter that converts a caloCand to a gctCand
void convertToGct(EmInputCandVec cand);

int main()
{
   EmInputCandVec inputs;
  EmInputCandVec sourceIso_0; 
  EmInputCandVec sourceIso_1;
  EmInputCandVec sourceIso_2;
  EmOutputCandVec tempVec;
  EmOutputCandVec outputs;
 
  // Number of bunch crossings to look at
  int bx = 1;

  //Construct 3 source cards (1 Rct crate)
  L1GctSourceCard* source_0 = new L1GctSourceCard(0,L1GctSourceCard::cardType1);
  L1GctSourceCard* source_1 = new L1GctSourceCard(1,L1GctSourceCard::cardType1);
  L1GctSourceCard* source_2 = new L1GctSourceCard(2,L1GctSourceCard::cardType1);

  //Constructor with 6 non iso electron candidates
  L1GctElectronSorter* testSort = new L1GctElectronSorter(4*bx,1);
  const string testFile_0 = "data/testElectronsRct_0";
  const string testFile_1 = "data/testElectronsRct_1";
  const string testFile_2 = "data/testElectronsRct_2";
 
  source_0->openInputFile(testFile_0);
  source_1->openInputFile(testFile_1);
  source_2->openInputFile(testFile_2);
  source_0->readBX();
  source_1->readBX();
  source_2->readBX();
  sourceIso_0 = source_0->getNonIsoElectrons();
  sourceIso_1 = source_1->getNonIsoElectrons(); 
  sourceIso_2 = source_2->getNonIsoElectrons();

  //Load in the data files for comparison
  LoadFileData(testFile_0,1);
  EmOutputCandVec fileData_0; 
  LoadFileData(testFile_1,1);
  EmOutputCandVec fileData_1; 
  LoadFileData(testFile_2,1);
  EmOutputCandVec fileData_2 = gctData;
  cout<<"\n ================= From data files ======================"<<endl;
  print(fileData_0);
  print(fileData_1);
  print(fileData_2);
  cout << "\n ================= From source cards ====================="<<endl;
  convertToGct(sourceIso_0);
  print(gctData);
  convertToGct(sourceIso_1);
  print(gctData);
  convertToGct(sourceIso_2);
  print(gctData);
  cout <<"\n ====== Eta phi table for electrons in the Rct crates ====="<<endl;
  eta_phiTable(1);
  eta_phiTable(10);

  //clean up
  delete testSort;
  delete source_0;
  delete source_1;
  delete source_2;
  return 0;   
}

// Function that outputs the eta and phi values for a given rct crate
void eta_phiTable(int crateNo){
  vector<float> phi(18);
  vector<float> eta(7);
  vector<int> card(6);

  //Fill vectors with appropiate values                                                                 
  for(int i=0;i<7;i++){
    eta[i] = i;
    card[i] = i;
  }
  
  for(int i=0;i<18;i++){
    phi[i] = i;
  }
  
  //Print values                                                                                        
  for(int i=0;i<18;i++){ //Loop over crates                                                             
    if(i==crateNo){
    cout <<"======================================================================" <<endl;
    cout <<"For crate["<<std::dec<<i<<"] "<<endl;
    
      for(int n=0;n<7;n++){ //Loop over cards                                                             
	int signEta = 1;
	int crate = 18;
	if(i<9){
	  signEta = -1;
	  crate = 0;
	}
	for(int m=0;m<2;m++){ //Loop over regions                                                         
	  if(n<3 || (n==6 && m==0)){
	    cout <<" Card["<<n<<"]: Region["<<m<<"] Phi is "<<phi[i]+1+i-crate<<" and Eta is ";
	    if(n==6){
	      cout<<eta[n]*signEta<<endl;
	    }else{
	      cout<<(eta[n]+m+n)*signEta<<endl;
	    }
	  }else{
	    cout <<" Card["<<n<<"]: Region["<<m<<"] Phi is "<<phi[i]+i-crate<<" and Eta is ";
	    if(n==6 && m==1){
	      cout<<eta[n]*signEta<<endl;
	    }else{
	      cout<<(eta[n-3]+m+(n-3))*signEta<<endl;
	    }
	  }
	}
      }
    }
  }
  return;
}

void LoadFileData(const string &inputFile, bool iso)
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
  //for(int n=0; n<bxs; n++){ //Loop over the bx bunch crossings                  
    file>> bxNo;
    file >>std::dec>> dummy;
    for(int i=0; i<58; i++){ //Loops over one bunch-crossing                    
      if(i<8){
        if(i>3 && iso){
          file >>std::hex>> dummy;
          candRank = dummy & 0x3f;
          candRegion = (dummy>>6) & 0x1;
          candCard = (dummy>>7) & 0x7;
          candIso = 1;
          L1CaloEmCand electrons(candRank, candRegion, candCard, candCrate, candIso);
	  data.push_back(electrons);
        }else{
          if(i<4 && !iso){
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
  convertToGct(data);
  file.close();
  return;
}

void print(EmOutputCandVec candidates){
  EmOutputCandVec cands = candidates; 
  for(unsigned int i=0; i!=cands.size(); i++){
    cout<<"          Rank: "<<cands[i].rank()<<"  Eta: "<<cands[i].etaIndex()<<"  Phi: "<<cands[i].phiIndex()<<endl;
  }
  return;
}

void convertToGct(vector<L1CaloEmCand> cand)
{
  gctData.resize(cand.size());
  for(unsigned int i = 0;i!=cand.size();i++){
    unsigned rank = cand[i].rank();
    unsigned card = cand[i].rctCard();
    unsigned region = cand[i].rctRegion();
    //unsigned crate = cand[i].rctCrate();
    //bool sign = (crate<9?1:0);
    bool isolation = cand[i].isolated();
    unsigned eta = 10; //initialisation value, outside eta range                     
    unsigned phi = 50;

    //    cout<<"Crate["<<crate<<"] with rank: "<<rank<<", card: "<<card<<" and region: "<<region<<" and sign: "<<sign<<endl;

    switch(card){
    case 0:
      phi = 1;
      if(region == 0){
        eta = 0;
      }else{
        eta = 1;
      }
      //      cout<<"Switch 0 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    case 1:
      phi = 1;
      if(region == 0){
        //cout<<"Switch 1 for phi: "<<phi<<" and eta: "<<eta<<endl;
        eta = 2;
      }else{
        eta = 3;
      }
      //cout<<"Switch 1 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    case 2:
      phi = 1;
      if(region == 0){
        eta = 4;
      }else{
        eta = 5;
      }
      //cout<<"Switch 2 for phi: "<<phi<<" and eta: "<<eta<<endl;
    case 3:
      phi = 0;
      if(region == 0){
        eta = 0;
      }else{
        eta = 1;
      }
      //cout<<"Switch 3 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    case 4:
      phi = 0;
      if(region == 0){
        eta = 2;
      }else{
        eta = 3;
      }
      //      cout<<"Switch 4 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    case 5:
      phi = 0;
      if(region == 0){
        eta = 4;
      }else{
        eta = 5;
      }
      //cout<<"Switch 5 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    case 6:
      if(region == 0){
        eta = 6;
        phi = 1;
      }else{
        eta = 6;
        phi = 0;
      }
      //cout<<"Switch 6 for phi: "<<phi<<" and eta: "<<eta<<endl;
      break;
    }
    L1GctEmCand gctTemp(rank,phi,eta,isolation);
    gctData[i] = gctTemp;

  }

  return;
}
