#include "L1Trigger/GlobalCaloTrigger/test/L1GctCompare.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1GctCompare::L1GctCompare(const edm::ParameterSet& iConfig) :
   m_cJets_tag1(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi1","l1GctEmulDigis"),"cenJets"),
   m_cJets_tag2(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi2","l1GctEmulDigi2"),"cenJets"),
   m_tJets_tag1(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi1","l1GctEmulDigis"),"tauJets"),
   m_tJets_tag2(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi2","l1GctEmulDigi2"),"tauJets"),
   m_fJets_tag1(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi1","l1GctEmulDigis"),"forJets"),
   m_fJets_tag2(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi2","l1GctEmulDigi2"),"forJets"),
   m_energy_tag1(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi1","l1GctEmulDigis")),
   m_energy_tag2(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigi2","l1GctEmulDigi2")),
   m_jets(0.0), m_totDiff(0.0), m_totDiffSq(0.0)
{
   //now do what ever initialization is needed


}


L1GctCompare::~L1GctCompare()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1GctCompare::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  
  // Get the L1 candidates from the event
  Handle<L1GctJetCandCollection> cJets1;
  iEvent.getByLabel(m_cJets_tag1,cJets1);

  Handle<L1GctJetCandCollection> cJets2;
  iEvent.getByLabel(m_cJets_tag2,cJets2);

  Handle<L1GctJetCandCollection> tJets1;
  iEvent.getByLabel(m_tJets_tag1,tJets1);

  Handle<L1GctJetCandCollection> tJets2;
  iEvent.getByLabel(m_tJets_tag2,tJets2);

  Handle<L1GctJetCandCollection> fJets1;
  iEvent.getByLabel(m_fJets_tag1,fJets1);

  Handle<L1GctJetCandCollection> fJets2;
  iEvent.getByLabel(m_fJets_tag2,fJets2);

  if (cJets1->size()==4 && cJets2->size()==4) { 
    for (size_t i=0; i<4; i++) {
      if (((*cJets1).at(i).etaIndex() == (*cJets2).at(i).etaIndex()) &&
          ((*cJets1).at(i).etaSign()  == (*cJets2).at(i).etaSign()) &&
          ((*cJets1).at(i).phiIndex() == (*cJets2).at(i).phiIndex())) {
        if (((*cJets1).at(i).rank()>0) || ((*cJets2).at(i).rank()>0)) {
          float r1 = static_cast<float>((*cJets1).at(i).rank());
          float r2 = static_cast<float>((*cJets2).at(i).rank());
          theJetRankComparison->Fill(r1,r2);
          theJetRankDifference->Fill((r1-r2),r1);
          theCenJetRankComparison->Fill(r1,r2);
          theCenJetRankDifference->Fill((r1-r2),r1);
          unsigned eta = (*cJets1).at(i).regionId().rctEta();
          theJetRankDifferencePerEtaBin.at(eta)->Fill((r1-r2),r1);
          m_jets      += 1.0;
          m_totDiff   += (r1-r2);
          m_totDiffSq += (r1-r2)*(r1-r2);
        }
      }
    }
  }

  if (tJets1->size()==4 && tJets2->size()==4) {
    for (size_t i=0; i<4; i++) {
      if (((*tJets1).at(i).etaIndex() == (*tJets2).at(i).etaIndex()) &&
          ((*tJets1).at(i).etaSign()  == (*tJets2).at(i).etaSign()) &&
          ((*tJets1).at(i).phiIndex() == (*tJets2).at(i).phiIndex())) {
        if (((*tJets1).at(i).rank()>0) || ((*tJets2).at(i).rank()>0)) {
          float r1 = static_cast<float>((*tJets1).at(i).rank());
          float r2 = static_cast<float>((*tJets2).at(i).rank());
          theJetRankComparison->Fill(r1,r2);
          theJetRankDifference->Fill((r1-r2),r1);
          theTauJetRankComparison->Fill(r1,r2);
          theTauJetRankDifference->Fill((r1-r2),r1);
          unsigned eta = (*tJets1).at(i).regionId().rctEta();
          theJetRankDifferencePerEtaBin.at(eta)->Fill((r1-r2),r1);
          m_jets      += 1.0;
          m_totDiff   += (r1-r2);
          m_totDiffSq += (r1-r2)*(r1-r2);
        }
      }
    }
  }
  if (fJets1->size()==4 && fJets2->size()==4) {
    for (size_t i=0; i<4; i++) {
      if (((*fJets1).at(i).etaIndex() == (*fJets2).at(i).etaIndex()) &&
          ((*fJets1).at(i).etaSign()  == (*fJets2).at(i).etaSign()) &&
          ((*fJets1).at(i).phiIndex() == (*fJets2).at(i).phiIndex())) {
        if (((*fJets1).at(i).rank()>0) || ((*fJets2).at(i).rank()>0)) {
          float r1 = static_cast<float>((*fJets1).at(i).rank());
          float r2 = static_cast<float>((*fJets2).at(i).rank());
          theJetRankComparison->Fill(r1,r2);
          theJetRankDifference->Fill((r1-r2),r1);
          theFwdJetRankComparison->Fill(r1,r2);
          theFwdJetRankDifference->Fill((r1-r2),r1);
          unsigned eta = (*fJets1).at(i).regionId().rctEta();
          theJetRankDifferencePerEtaBin.at(eta)->Fill((r1-r2),r1);
          m_jets      += 1.0;
          m_totDiff   += (r1-r2);
          m_totDiffSq += (r1-r2)*(r1-r2);
        }
      }
    }
  }

  // get the L1 energy sums from the event
  Handle< L1GctEtTotal > sumEt1 ;
  iEvent.getByLabel( m_energy_tag1, sumEt1 ) ;
  Handle< L1GctEtTotal > sumEt2 ;
  iEvent.getByLabel( m_energy_tag2, sumEt2 ) ;

  Handle< L1GctEtHad > sumHt1 ;
  iEvent.getByLabel( m_energy_tag1, sumHt1 ) ;
  Handle< L1GctEtHad > sumHt2 ;
  iEvent.getByLabel( m_energy_tag2, sumHt2 ) ;

  Handle< L1GctEtMiss > missEt1 ;
  iEvent.getByLabel( m_energy_tag1, missEt1 ) ;
  Handle< L1GctEtMiss > missEt2 ;
  iEvent.getByLabel( m_energy_tag2, missEt2 ) ;

  float etTot1  = static_cast<float>(sumEt1->et());
  float etTot2  = static_cast<float>(sumEt2->et());
  float etHad1  = static_cast<float>(sumHt1->et());
  float etHad2  = static_cast<float>(sumHt2->et());
  float etMiss1 = static_cast<float>(missEt1->et());
  float etMiss2 = static_cast<float>(missEt2->et());
  float etMPhi1 = static_cast<float>(missEt1->phi());
  float etMPhi2 = static_cast<float>(missEt2->phi());

  theSumEtComparison->Fill ( etTot1,  etTot2 );
  theSumEtDifference->Fill ( etTot1 - etTot2 );
  theSumHtComparison->Fill ( etHad1,  etHad2 );
  theSumHtDifference->Fill ( etHad1 - etHad2 );
  theMissEtComparison->Fill(etMiss1,  etMiss2);
  theMissEtDifference->Fill(etMiss1 - etMiss2);
  theMissEtPhiGCT1->Fill(etMPhi1);
  theMissEtPhiGCT2->Fill(etMPhi2);

  theSumHtRatio->Fill(etHad1/etHad2);

}


// ------------ method called once each job just before starting event loop  ------------
void 
L1GctCompare::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  // Comparison histograms for jet rank
  TFileDirectory dir = fs->mkdir("L1JetCompare");

  theJetRankComparison    = dir.make<TH2F>("JetRankComparison",    "Jet Rank gct1 vs gct2",
                                          64, 0., 64., 64, 0., 64.);
  theJetRankDifference    = dir.make<TH2F>("JetRankDifference",    "Jet Rank gct1-gct2 vs gct1",
                                          40, -20., 20., 64, 0., 64.);
  theCenJetRankComparison = dir.make<TH2F>("CenJetRankComparison", "Jet Rank gct1 vs gct2",
                                          64, 0., 64., 64, 0., 64.);
  theCenJetRankDifference = dir.make<TH2F>("CenJetRankDifference", "Jet Rank gct1-gct2 vs gct1",
                                          40, -20., 20., 64, 0., 64.);
  theTauJetRankComparison = dir.make<TH2F>("TauJetRankComparison", "Jet Rank gct1 vs gct2",
                                          64, 0., 64., 64, 0., 64.);
  theTauJetRankDifference = dir.make<TH2F>("TauJetRankDifference", "Jet Rank gct1-gct2 vs gct1",
                                          40, -20., 20., 64, 0., 64.);
  theFwdJetRankComparison = dir.make<TH2F>("FwdJetRankComparison", "Jet Rank gct1 vs gct2",
                                          64, 0., 64., 64, 0., 64.);
  theFwdJetRankDifference = dir.make<TH2F>("FwdJetRankDifference", "Jet Rank gct1-gct2 vs gct1",
                                          40, -20., 20., 64, 0., 64.);

  // Comparison histograms for jet rank for each eta bin
  TFileDirectory dir1 = fs->mkdir("L1JetsPerEtaBin");

  for (unsigned eta=0; eta<11; eta++) {
    std::stringstream ss;
    std::string title;
    ss << "JetRankDifferenceEtaBin" << eta;
    ss >> title;
    theJetRankDifferencePerEtaBin.push_back(dir1.make<TH2F>(title.c_str(), "Jet Rank gct1-gct2 vs gct1",
                                                            40, -20., 20., 64, 0., 64.));
  }

  // Comparison histograms for energy sums and missing Et
  TFileDirectory dir2 = fs->mkdir("L1EnergySums");
  theSumEtComparison  = dir2.make<TH2F>("SumEtComparison",    "Total Et gct1 vs gct2",
                                       128, 0., 512., 128, 0., 512.);
  theSumHtComparison  = dir2.make<TH2F>("SumHtComparison",    "Total Ht gct1 vs gct2",
                                       128, 0.,1024., 128, 0., 512.);
  theMissEtComparison = dir2.make<TH2F>("MissEtComparison", "Missing Et gct1 vs gct2",
                                       128, 0., 256., 128, 0., 256.);

  theSumEtDifference  = dir2.make<TH1F>("SumEtDifference",    "Total Et gct1-gct2", 40, -20., 20.);
  theSumHtDifference  = dir2.make<TH1F>("SumHtDifference",    "Total Ht gct1-gct2", 40, -20., 20.);
  theMissEtDifference = dir2.make<TH1F>("MissEtDifference", "Missing Et gct1-gct2", 40, -20., 20.);

  theMissEtPhiGCT1 = dir2.make<TH1F>("MissEtPhiGCT1", "Missing Et angle gct1", 72, 0., 72.);
  theMissEtPhiGCT2 = dir2.make<TH1F>("MissEtPhiGCT2", "Missing Et angle gct2", 72, 0., 72.);

  theSumHtRatio       = dir2.make<TH1F>("SumHtRatio",         "Total Ht gct1/gct2", 200, 0., 4.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1GctCompare::endJob() {
  float aveDiff = m_totDiff/m_jets;
  float varDiff = m_totDiffSq/m_jets - aveDiff*aveDiff;
  edm::LogInfo("L1GctCompare") << "Mean difference between jet ranks " << aveDiff
                               << " variance " << varDiff << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1GctCompare);
