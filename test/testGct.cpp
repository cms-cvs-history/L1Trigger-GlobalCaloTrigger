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
  cout<<"In main"<<endl;
  L1GlobalCaloTrigger* gct = L1GlobalCaloTrigger::theGct();
  //vector<L1GctEmCand> electrons = gct->getIsoElectrons();
  // vector<L1GctSourceCard*> theSourceCards;
  //theSourceCards = gct->getSourceCards();
  cout<<"Works so far"<<endl;
  //for(unsigned int i=0;i!=theSourceCards.size();i++){
    ///    theSourceCards[i]->openInputFile(dummyData);
    //vector<L1GctEmCand> isoElectrons = theSourceCards[i]->getIsoElectrons();
    //vector<L1GctEmCand> nonIsoElectrons = theSourceCards[i]->getNonIsoElectrons();
    //cout<<"For source card no: "<<i;
    //for(unsigned int n=0;n!=isoElectrons.size();n++){
  //cout<<" iso electrons are:"<<endl;
  //   cout<<isoElectrons[n].getRank()<<" "<<isoElectrons[n].getEta()<<" "<<isoElectrons[n].getPhi()<<endl;
  //  cout<<"and non-iso electrons are: "<<endl;
  //   cout<<nonIsoElectrons[n].getRank()<<" "<<isoElectrons[n].getEta()<<" "<<isoElectrons[n].getPhi()<<endl;
  //}
    //  Gct->process();
  //}
}
