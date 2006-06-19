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
 
    // Open source card input files
    cout << "Opening input source card files" << endl;
    gct->openSourceCardFiles("data/testElectronsRct_"); 

    // clear everything
    cout << "Reset the GCT" << endl; 
    gct->reset(); 

    // Run
    cout << "Run the GCT...." << endl;
    gct->process(); // process

    // central jet outputs to GT
    vector<L1GctJet> centralJets = gct->getCentralJets();
  
    // forward jet outputs to GT
    vector<L1GctJet> forwardJets = gct->getForwardJets();
  
    // tau jet outputs to GT
    vector<L1GctJet> tauJets = gct->getTauJets();

    // Now take a look at the output jets

    cout << "===== CENTRAL JETS ====" << endl;
    for (unsigned i=0; i<centralJets.size();i++){
      cout << centralJets.at(i);
    }

    cout << "===== FORWARD JETS ====" << endl;
    for (unsigned i=0; i<forwardJets.size();i++){
      cout << forwardJets.at(i);
    }

    cout << "===== TAU JETS ====" << endl;
    for (unsigned i=0; i<tauJets.size();i++){
      cout << tauJets.at(i);
    }
    
  }
  catch (cms::Exception& e)
    {
      cerr << "CMS exception from " << e.what() << endl;
    }
  catch (std::exception& e)
    {
      cerr << "std exception from " << e.what() << endl;
    }
  catch (...) { // Catch anything unknown that goes wrong! 
    cerr << "yikes! something awful happened.... (unknown exception)" << endl; 
  } 
}

