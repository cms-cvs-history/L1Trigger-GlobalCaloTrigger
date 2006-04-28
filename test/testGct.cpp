/*! \file testGct.cpp
 * \test file for testing the (eventually) the full Gct 
 *
 *  For now, this test program only handles electrons.
 *  The 54 source cards handles 3 RCT crates each. The first RCT crate on
 *  each source card carries the electrons.
 *
 * \author Maria Hansen
 * \date April 2006
 */


#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEmCand.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

using namespace std;

int main()
{
 
  L1GlobalCaloTrigger* gct = L1GlobalCaloTrigger::theGct();
  //Firstly, check that there's nothing in the buffers of the GCT
  vector<L1GctEmCand> electrons = gct->getIsoElectrons();
  vector<L1GctEmCand> nonIso = gct->getNonIsoElectrons();
  cout<<"From GCT: iso electrons     &    non-iso electrons:"<<endl;
  cout<<"          Rank   Eta   Phi       Rank   Eta   Phi" <<endl;
  for(unsigned int i=0;i!=electrons.size();i++){
    cout<<"          "<<electrons[i].getRank()<<"      "<<electrons[i].getEta()<<"     "<<electrons[i].getPhi()<<
    "         "<<nonIso[i].getRank()<<"      "<<nonIso[i].getEta()<<"     "<<nonIso[i].getPhi()<<endl;
  }

  //Open a source card file and look at the electrons in there
  vector<L1GctSourceCard*> theSourceCards = gct->getSourceCards();
  for(unsigned int i=0;i!=(theSourceCards.size()/3);i++){
    //These ready for later when gotten more than one BX/file to read in
    //  theSourceCards[i]->openInputFile("testSourceCardInput.txt");
    //  theSourceCards[i*3]->readBX();
    //First checking buffers before reading in data
  vector<L1GctEmCand> isoElectrons = theSourceCards[i*3]->getIsoElectrons();
  vector<L1GctEmCand> nonIsoElectrons = theSourceCards[i*3]->getNonIsoElectrons();
  cout<<"For source card no: "<<i<<"  Rank   Eta   Phi"<<endl;
  for(unsigned int n=0;n!=isoElectrons.size();n++){
    cout<<"Iso electrons are:      "<<isoElectrons[n].getRank()<<"      "<<isoElectrons[n].getEta()<<"     "<<isoElectrons[n].getPhi()<<endl;
    cout<<"Non-iso are:            "<<nonIsoElectrons[n].getRank()<<"      "<<isoElectrons[n].getEta()<<"     "<<isoElectrons[n].getPhi()<<endl;
  }
  }
  theSourceCards[0]->openInputFile("testSourceCardInput.txt");
  theSourceCards[0]->readBX();
  vector<L1GctEmCand> iso_zero = theSourceCards[0]->getIsoElectrons();
  vector<L1GctEmCand> nonIso_zero = theSourceCards[0]->getNonIsoElectrons();
  cout<<"For source card no  0   Rank   Eta   Phi"<<endl;
   for(unsigned int n=0;n!=iso_zero.size();n++){
     cout<<"Iso electrons are:      "<<std::hex<<iso_zero[n].getRank()<<"      "<<std::hex<<iso_zero[n].getEta()<<"     "<<std::hex<<iso_zero[n].getPhi()<<endl;
      cout<<"Non-iso are:            "<<nonIso_zero[n].getRank()<<"      "<<nonIso_zero[n].getEta()<<"     "<<nonIso_zero[n].getPhi()<<endl;
   }
}

