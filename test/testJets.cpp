/*! \file testJets.cpp
 * \test file for testing the full GCT class jet algorithms
 *
 *  
 *
 * \author Alex Tapper
 * \date May 2006
 */

#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h" 
#include "FWCore/Utilities/interface/Exception.h"

#include <iostream>
#include <exception>
#include <vector>

using namespace std;

int main()
{
  try {
    // New GCT
    L1GlobalCaloTrigger* gct = new L1GlobalCaloTrigger(true); 
    cout << "Reset the GCT" << endl;

    // clear everything
    gct->reset(); 

    // Open source card input files
    cout << "Opening input source card files" << endl;
    gct->openSourceCardFiles("RCT_"); 

    // Run
    cout << "Run the GCT...." << endl;
    gct->process(); // process

    // central jet outputs to GT
    vector<L1GctJetCand> centralJets = gct->getCentralJets();
  
    // forward jet outputs to GT
    vector<L1GctJetCand> forwardJets = gct->getForwardJets();
  
    // tau jet outputs to GT
    vector<L1GctJetCand> tauJets = gct->getTauJets();

    // Now take a look at the output jets

    cout << "===== CENTRAL JETS ====" << endl;
    for (unsigned i=0; i<centralJets.size();i++){
      cout << centralJets[i];
    }

    cout << "===== FORWARD JETS ====" << endl;
    for (unsigned i=0; i<forwardJets.size();i++){
      cout << forwardJets[i];
    }

    cout << "===== TAU JETS ====" << endl;
    for (unsigned i=0; i<tauJets.size();i++){
      cout << tauJets[i];
    }
    
  }
  catch (cms::Exception& e)
    {
      cerr << e.what() << endl;
    }
  catch (...) { // Catch anything unknown that goes wrong! 
    cerr << "yikes! something awful happened.... (unknown exception)" << endl; 
  } 
}

