#ifndef L1GCTVALIDATION_H
#define L1GCTVALIDATION_H
// -*- C++ -*-
//
// Package:    L1GlobalCaloTrigger
// Class:      L1GctValidation
// 
/**\class L1GctValidation L1GctValidation.cc L1Trigger/L1GlobalCaloTrigger/plugins/L1GctValidation.cc

 Description: produces standard plots of Gct output quantities to enable validation
              of global event quantities in particular

*/
//
// Author: Greg Heath
// Date:   February 2008
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "TH1.h"
#include "TH2.h"
//
// class declaration
//

class L1GctValidation : public edm::EDAnalyzer {
   public:
      explicit L1GctValidation(const edm::ParameterSet&);
      ~L1GctValidation();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      edm::InputTag m_energy_tag;
      edm::InputTag m_missHt_tag;

      TH1F* theSumEtInLsb;
      TH1F* theSumHtInLsb;
      TH1F* theMissEtInLsb;
      TH1F* theMissHtInLsb;
      TH1F* theSumEtInGeV;
      TH1F* theSumHtInGeV;
      TH1F* theMissEtInGeV;
      TH1F* theMissEtAngle;
      TH2F* theMissEtVector;
      TH1F* theMissHtInGeV;
      TH1F* theMissHtAngle;
      TH2F* theMissHtVector;

      TH2F* theMissEtVsMissHt;
      TH2F* theMissEtVsMissHtAngle;

      std::vector<TH1F*> theJetCounts;
      TH1F* theHfRing0EtSumPositiveEta;
      TH1F* theHfRing0EtSumNegativeEta;
      TH1F* theHfRing1EtSumPositiveEta;
      TH1F* theHfRing1EtSumNegativeEta;
      TH1F* theHfRing0CountPositiveEta;
      TH1F* theHfRing0CountNegativeEta;
      TH1F* theHfRing1CountPositiveEta;
      TH1F* theHfRing1CountNegativeEta;

};
#endif