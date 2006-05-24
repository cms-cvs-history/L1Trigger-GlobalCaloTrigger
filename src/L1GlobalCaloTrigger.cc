#include "L1Trigger/GlobalCaloTrigger/interface/L1GlobalCaloTrigger.h"

#include "L1Trigger/GlobalCaloTrigger/interface/L1GctSourceCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetLeafCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctEmLeafCard.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelJetFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctWheelEnergyFpga.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinalStage.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctGlobalEnergyAlgos.h"
#include "L1Trigger/GlobalCaloTrigger/interface/L1GctElectronFinalSort.h"

#include <sstream>

using std::cout;
using std::string;
using std::stringstream;
using std::endl;
using std::vector;

// constructor
L1GlobalCaloTrigger::L1GlobalCaloTrigger() :
  theSourceCards(54),
  theJetLeafCards(6),
  theEmLeafCards(2),
  theWheelJetFpgas(2),
  theWheelEnergyFpgas(2)
{
  
  build();
  setup();
  
}

L1GlobalCaloTrigger::~L1GlobalCaloTrigger()
{
  theSourceCards.clear();
}

void L1GlobalCaloTrigger::openSourceCardFiles(string fileBase){
 //Loop running over the 18 RCT-crate files, allocating 3 sourcecards per file
 for(int i = 0;i < 18; i++){
   string fileNo;
   stringstream ss;
   ss << i;
   ss >> fileNo;
   string fileName = fileBase+fileNo;
   theSourceCards[3*i]->openInputFile(fileName);
   theSourceCards[3*i+1]->openInputFile(fileName);
   theSourceCards[3*i+2]->openInputFile(fileName);
 }
}

void L1GlobalCaloTrigger::reset() {

  // Source cards
  for (int i=0; i<54; i++) {
    theSourceCards[i]->reset();
  }

  // EM Leaf Card
  for (int i=0; i<4; i++) {
    theEmLeafCards[i]->reset();
  }

  // Jet Leaf cards
  for (int i=0; i<6; i++) {
    theJetLeafCards[i]->reset();
  }

  // Wheel Cards
  for (int i=0; i<2; i++) {
    theWheelJetFpgas[i]->reset();
  }

  for (int i=0; i<2; i++) {
    theWheelEnergyFpgas[i]->reset();
  }

  // Electron Final Stage
  theIsoEmFinalStage->reset();
  theNonIsoEmFinalStage->reset();

  // Jet Final Stage
  theJetFinalStage->reset();

  // Energy Final Stage
  theEnergyFinalStage->reset();

}

void L1GlobalCaloTrigger::process() {
		
  // Source cards
  for (int i=0; i<54; i++) {
    theSourceCards[i]->fetchInput();
  }

  // EM Leaf Card
  for (int i=0; i<4; i++) {
    theEmLeafCards[i]->fetchInput();
    theEmLeafCards[i]->process();
  }

  // Jet Leaf cards
  for (int i=0; i<6; i++) {
    theJetLeafCards[i]->fetchInput();
    theJetLeafCards[i]->process();
  }

  // Wheel Cards
  for (int i=0; i<2; i++) {
    theWheelJetFpgas[i]->fetchInput();
    theWheelJetFpgas[i]->process();
  }

  for (int i=0; i<2; i++) {
    theWheelEnergyFpgas[i]->fetchInput();
    theWheelEnergyFpgas[i]->process();
  }

  // Electron Final Stage
  theIsoEmFinalStage->fetchInput();
  theIsoEmFinalStage->process();

  theNonIsoEmFinalStage->fetchInput();
  theNonIsoEmFinalStage->process();


  // Jet Final Stage
  theJetFinalStage->fetchInput();
  theJetFinalStage->process();

  // Energy Final Stage
  theEnergyFinalStage->fetchInput();
  theEnergyFinalStage->process();
	
}

void L1GlobalCaloTrigger::print() {

	cout << "=== Global Calo Trigger ===" << endl;
	cout << "=== START DEBUG OUTPUT  ===" << endl;

	cout << endl;
	cout << "N Source Cards " << theSourceCards.size() << endl;
	cout << "N Jet Leaf Cards " << theJetLeafCards.size() << endl;
	cout << "N Wheel Jet Fpgas " << theWheelJetFpgas.size() << endl;
	cout << "N Wheel Energy Fpgas " << theWheelEnergyFpgas.size() << endl;
	cout << "N Em Leaf Cards" << theEmLeafCards.size() << endl;
	cout << endl;

	for (unsigned i=0; i<theSourceCards.size(); i++) {
          cout << "Source Card : " << i <<endl;
          cout << (*theSourceCards[i]); 
	}
	cout << endl;

	for (unsigned i=0; i<theJetLeafCards.size(); i++) {
	  cout << "Jet Leaf : " << i << endl;
          cout << (*theJetLeafCards[i]);
	}
	cout << endl;

 	for (unsigned i=0; i<theWheelJetFpgas.size(); i++) {
          cout << "Wheel Jet FPGA : " << i << endl; 
          cout << (*theWheelJetFpgas[i]);
 	}
	cout << endl;

 	for (unsigned i=0; i<theWheelEnergyFpgas.size(); i++) {
          cout << "Wheel Energy FPGA : " << i << endl; 
	  cout << (*theWheelEnergyFpgas[i]);
 	}
	cout << endl;

 	cout << (*theJetFinalStage);
	cout << endl;

 	cout << (*theEnergyFinalStage);
	cout << endl;

	for (unsigned i=0; i<theEmLeafCards.size(); i++) {
	  cout << "EM Leaf : " << i << endl;
          cout << (*theEmLeafCards[i]);
	}
	cout << endl;

 	cout << (*theIsoEmFinalStage);
	cout << endl;

	cout << (*theNonIsoEmFinalStage);

	cout << "=== Global Calo Trigger ===" << endl;
	cout << "===  END DEBUG OUTPUT   ===" << endl;
 
}

// isolated EM outputs
vector<L1GctEmCand> L1GlobalCaloTrigger::getIsoElectrons() { 
  return theIsoEmFinalStage->OutputCands();
}	

// non isolated EM outputs
vector<L1GctEmCand> L1GlobalCaloTrigger::getNonIsoElectrons() {
return theNonIsoEmFinalStage->OutputCands(); 
}

// central jet outputs to GT
vector<L1GctJetCand> L1GlobalCaloTrigger::getCentralJets() {
 return theJetFinalStage->getCentralJets();
}

// forward jet outputs to GT
vector<L1GctJetCand> L1GlobalCaloTrigger::getForwardJets() { 
  return theJetFinalStage->getForwardJets(); 
}

// tau jet outputs to GT
vector<L1GctJetCand> L1GlobalCaloTrigger::getTauJets() { 
  return theJetFinalStage->getTauJets(); 
}

// total Et output
unsigned L1GlobalCaloTrigger::getEtSum() {
  return 0; //theEnergyFinalStage->getEtSum();
}

unsigned L1GlobalCaloTrigger::getEtHad() {
  return 0; //theEnergyFinalStage->getEtHad();
}

unsigned L1GlobalCaloTrigger::getEtMiss() {
  return 0; //theEnergyFinalStage->getEtMiss();
}

unsigned L1GlobalCaloTrigger::getEtMissPhi() {
  return 0; //theEnergyFinalStage->etMissPhi();
}





/* PRIVATE METHODS */

// instantiate hardware/algorithms
void L1GlobalCaloTrigger::build() {

  for (int i=0; i<18; i++) {
    theSourceCards[3*i] = new L1GctSourceCard(3*i, L1GctSourceCard::cardType1);
    theSourceCards[3*i+1] = new L1GctSourceCard(3*i+1, L1GctSourceCard::cardType2);
    theSourceCards[3*i+2] = new L1GctSourceCard(3*i+2, L1GctSourceCard::cardType3);
  }

  for (int i=0; i<6; i++) {
    theJetLeafCards[i] = new L1GctJetLeafCard(i, i % 3);
  }
  
  for (int i=0; i<2; i++) {
    theEmLeafCards[i] = new L1GctEmLeafCard(i);
  }
    
  for (int i=0; i<2; i++) {
    theWheelJetFpgas[i] = new L1GctWheelJetFpga(i);
    theWheelEnergyFpgas[i] = new L1GctWheelEnergyFpga(i);
  }
  
  theJetFinalStage = new L1GctJetFinalStage();
  theEnergyFinalStage = new L1GctGlobalEnergyAlgos();

  theIsoEmFinalStage = new L1GctElectronFinalSort(true);
  theNonIsoEmFinalStage = new L1GctElectronFinalSort(false);

}

// wire up the hardware/algos
void L1GlobalCaloTrigger::setup() {

  // EM leaf cards
  for (int i=0; i<2; i++) {
    for (int j=0; j<9; j++) {
      theEmLeafCards[i]->setInputSourceCard(j, theSourceCards[(i*9+j)*3]);
    }
  }

  // Jet Leaf cards
  for (int i=0; i<6; i++) {
    for (int j=0; j<3; j++) {
      for (int k=1; k<3; k++) {
	theJetLeafCards[i]->setInputSourceCard(j, theSourceCards[(i*3+j)*3+k]);
      }
      // Neighbour connections
      int iup = (i*3+3) % 9;
      int idn = (i*3+8) % 9;
      int ii, i0, i1, i2, i3, i4, i5;
      if (i<3) {
	ii = iup;
	// Remaining connections for the TDR jetfinder only
	i0 = idn;
	i1 = idn+9;
	i2 = i*3+9;
	i3 = i*3+10;
	i4 = i*3+11;
	i5 = iup+9;
      } else {
	ii = iup+9;
	// Remaining connections for the TDR jetfinder only
	i0 = idn+9;
	i1 = idn;
	i2 = i*3-9;
	i3 = i*3-8;
	i4 = i*3-7;
	i5 = iup;
      }
      theJetLeafCards[i]->setInputSourceCard( 6, theSourceCards[ii*3+1]);
      theJetLeafCards[i]->setInputSourceCard( 7, theSourceCards[ii*3+2]);
      // Remaining connections for the TDR jetfinder only
      theJetLeafCards[i]->setInputSourceCard( 8, theSourceCards[i0*3+1]);
      theJetLeafCards[i]->setInputSourceCard( 9, theSourceCards[i0*3+2]);
      theJetLeafCards[i]->setInputSourceCard(10, theSourceCards[i1*3+1]);
      theJetLeafCards[i]->setInputSourceCard(11, theSourceCards[i2*3+1]);
      theJetLeafCards[i]->setInputSourceCard(12, theSourceCards[i3*3+1]);
      theJetLeafCards[i]->setInputSourceCard(13, theSourceCards[i4*3+1]);
      theJetLeafCards[i]->setInputSourceCard(14, theSourceCards[i5*3+1]);
      //


    }
  }
  
  // Wheel Fpgas
  for (int i=0; i<2; i++) {
    for (int j=0; j<3; j++) {
      theWheelJetFpgas[i]->setInputLeafCard(j, theJetLeafCards[i*2+j]);
      theWheelEnergyFpgas[i]->setInputLeafCard(j, theJetLeafCards[i*2+j]);
    }
  }

  // Electron Final Sort
  for (int i=0; i<2; i++) {
    theIsoEmFinalStage->setInputLeafCard(i, theEmLeafCards[i]);
    theNonIsoEmFinalStage->setInputLeafCard(i, theEmLeafCards[i]);
  }

  // Jet Final Stage
  for (int i=0; i<2; i++) {
    theJetFinalStage->setInputWheelJetFpga(i, theWheelJetFpgas[i]);
  }

  // Global Energy Algos
  theEnergyFinalStage->setMinusWheelEnergyFpga(theWheelEnergyFpgas[0]);
  theEnergyFinalStage->setPlusWheelEnergyFpga(theWheelEnergyFpgas[1]);
  theEnergyFinalStage->setMinusWheelJetFpga(theWheelJetFpgas[0]);
  theEnergyFinalStage->setPlusWheelJetFpga(theWheelJetFpgas[1]);

}
