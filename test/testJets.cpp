/*! \file testJets.cpp
 * \test file for testing the full GCT 
 *
 *  
 *
 * \author Alex Tapper
 * \date May 2006
 */


#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h" 
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEmCand.h"

//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>

int main()
{
 
  L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger();
  gct->print();
  //Firstly, check that there's nothing in the buffers of the GCT
//   std::vector<L1GctEmCand> electrons = gct->getIsoElectrons();
//   std::vector<L1GctEmCand> nonIso = gct->getNonIsoElectrons();
//   std::cout<<"From GCT: iso electrons     &    non-iso electrons:"<<std::endl;
//   std::cout<<"          Rank   Eta   Phi       Rank   Eta   Phi" <<std::endl;
//   for(unsigned int i=0;i!=electrons.size();i++){
//     std::cout<<"          "<<electrons[i].rank()<<"      "<<electrons[i].eta()<<"     "<<electrons[i].phi()<<
//     "         "<<nonIso[i].rank()<<"      "<<nonIso[i].eta()<<"     "<<nonIso[i].phi()<<std::endl;
//   }

//   //Open a source card file and look at the electrons in there
//   vector<L1GctSourceCard*> theSourceCards = gct->getSourceCards();
//   for(unsigned int i=0;i!=(theSourceCards.size()/3);i++){
//     //These ready for later when gotten more than one BX/file to read in
//     //  theSourceCards[i]->openInputFile("testSourceCardInput.txt");
//     //  theSourceCards[i*3]->readBX();
//     //First checking buffers before reading in data
//   vector<L1GctEmCand> isoElectrons = theSourceCards[i*3]->getIsoElectrons();
//   vector<L1GctEmCand> nonIsoElectrons = theSourceCards[i*3]->getNonIsoElectrons();
//   cout<<"For source card no: "<<i<<"  Rank   Eta   Phi"<<endl;
//   for(unsigned int n=0;n!=isoElectrons.size();n++){
//     cout<<"Iso electrons are:      "<<isoElectrons[n].rank()<<"      "<<isoElectrons[n].eta()<<"     "<<isoElectrons[n].phi()<<endl;
//     cout<<"Non-iso are:            "<<nonIsoElectrons[n].rank()<<"      "<<isoElectrons[n].eta()<<"     "<<isoElectrons[n].phi()<<endl;
//   }
//   }
//   theSourceCards[0]->openInputFile("testSourceCardInput.txt");
//   theSourceCards[0]->readBX();
//   vector<L1GctEmCand> iso_zero = theSourceCards[0]->getIsoElectrons();
//   vector<L1GctEmCand> nonIso_zero = theSourceCards[0]->getNonIsoElectrons();
//   cout<<"For source card no  0   Rank   Eta   Phi"<<endl;
//    for(unsigned int n=0;n!=iso_zero.size();n++){
//      cout<<"Iso electrons are:      "<<std::hex<<iso_zero[n].rank()<<"      "<<std::hex<<iso_zero[n].eta()<<"     "<<std::hex<<iso_zero[n].phi()<<endl;
//       cout<<"Non-iso are:            "<<nonIso_zero[n].rank()<<"      "<<nonIso_zero[n].eta()<<"     "<<nonIso_zero[n].phi()<<endl;
//    }
}

