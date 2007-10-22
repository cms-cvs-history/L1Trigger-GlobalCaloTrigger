// -*- C++ -*-
//
// Package:    L1GctCompare
// Class:      L1GctCompare
// 
/**\class L1GctCompare L1GctCompare.cc L1Trigger/L1GctCompare/src/L1GctCompare.cc

 Description: compares the GCT output after running twice with different setup info

*/
//
// Original Author:  Gregory Heath
//         Created:  Fri Jul 27 16:09:48 CEST 2007
// $Id: L1GctCompare.h,v 1.2 2007/10/04 10:17:28 heath Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"

#include "TH1.h"
#include "TH2.h"
//
// class declaration
//

class L1GctCompare : public edm::EDAnalyzer {
   public:
      explicit L1GctCompare(const edm::ParameterSet&);
      ~L1GctCompare();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

      void analyzeJets(const edm::Handle<L1GctJetCandCollection>& j1,
                       const edm::Handle<L1GctJetCandCollection>& j2);

      edm::InputTag m_rctInput_tag;
      edm::InputTag m_cJets_tag1;
      edm::InputTag m_cJets_tag2;
      edm::InputTag m_tJets_tag1;
      edm::InputTag m_tJets_tag2;
      edm::InputTag m_fJets_tag1;
      edm::InputTag m_fJets_tag2;
      edm::InputTag m_energy_tag1;
      edm::InputTag m_energy_tag2;

      TH2F* theJetRankComparison;
      TH2F* theJetRankDifference;
      TH2F* theCenJetRankComparison;
      TH2F* theCenJetRankDifference;
      TH2F* theTauJetRankComparison;
      TH2F* theTauJetRankDifference;
      TH2F* theFwdJetRankComparison;
      TH2F* theFwdJetRankDifference;

      TH2F* theEtaPhiMatchedJets;
      TH1F* theJetRankGCT1Only;
      TH1F* theJetRankGCT2Only;
      TH2F* theEtaPhiGCT1Only;
      TH2F* theEtaPhiGCT2Only;

      std::vector<TH2F*> theJetRankDifferencePerEtaBin;

      TH2F* theSumEtComparison;
      TH2F* theSumHtComparison;
      TH2F* theMissEtComparison;

      TH1F* theSumEtDifference;
      TH1F* theSumHtDifference;
      TH1F* theMissEtDifference;

      TH1F* theMissEtPhiGCT1;
      TH1F* theMissEtPhiGCT2;

      TH2F* theMissEtVecRCT;
      TH2F* theMissEtVecGCT1;
      TH2F* theMissEtVecGCT2;
      TH2F* theMETOffsetGCT1;
      TH2F* theMETOffsetGCT2;

      std::vector<TH1F*> theRctExPerEtaBin;
      std::vector<TH1F*> theRctEyPerEtaBin;

      TH1F* theSumHtRatio;

      float m_jets, m_totDiff, m_totDiffSq;

};

