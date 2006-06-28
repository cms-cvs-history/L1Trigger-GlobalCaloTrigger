/*! \file testElectrons.cpp
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


//Standard library headers
#include <fstream>   //for file IO
#include <string>
#include <vector>
#include <iostream>
#include <sstream>  //for int->char conversion

using std::cout;
using std::endl;

using namespace std;

int main()
{
  try { 
  L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(1);
  
  //gct->print();
  //Firstly, check that there's nothing in the buffers of the GCT
  vector<L1GctEmCand> electrons = gct->getIsoElectrons();
  vector<L1GctEmCand> nonIso = gct->getNonIsoElectrons();
  cout<<"From GCT: iso electrons     &    non-iso electrons:"<<endl;
  cout<<"          Rank   Eta   Phi       Rank   Eta   Phi" <<endl;
  for(unsigned int i=0;i!=electrons.size();i++){
    cout<<"          "<<electrons[i].rank()<<"      "<<electrons[i].level1EtaIndex()<<"     "<<electrons[i].level1PhiIndex()<<
    "         "<<nonIso[i].rank()<<"      "<<nonIso[i].level1EtaIndex()<<"     "<<nonIso[i].level1PhiIndex()<<endl;
  }

  //Open source card files and look at the electrons in there
  cout<<"testelectrons open file?"<<endl;
    std::string fileName = "data/testElectronsRct_";
    gct->openSourceCardFiles(fileName);
    cout<<"about to process..."<<endl;
    gct->process();
    //gct->print();
    vector<L1GctEmCand> newIso = gct->getIsoElectrons();
    vector<L1GctEmCand> newnonIso = gct->getNonIsoElectrons();
    cout<<"From GCT: iso electrons     &    non-iso electrons:"<<endl;
    cout<<"          Rank   Eta   Phi       Rank   Eta   Phi" <<endl;
    for(unsigned int i=0;i!=newIso.size();i++){
      cout<<"          "<<newIso[i].rank()<<"      "<<newIso[i].level1EtaIndex()<<"     "<<newIso[i].level1PhiIndex()<<
	"         "<<newnonIso[i].rank()<<"      "<<newnonIso[i].level1EtaIndex()<<"     "<<newnonIso[i].level1PhiIndex()<<endl;
    }
    cout<<"END OF TEST PROGRAM"<<endl;
  }
  catch(std::exception e){
    cerr << e.what() << endl;
  }
}

