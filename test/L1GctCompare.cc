#include "L1Trigger/GlobalCaloTrigger/test/L1GctCompare.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1GlobalCaloTrigger/interface/L1GctEtSums.h"

#include <math.h>

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
   m_rctInput_tag(iConfig.getUntrackedParameter<std::string>("L1RctEmulDigis","l1RctEmulDigis")),
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

  // Get the Rct input from the event
  Handle<L1CaloEmCollection> inputEmCands;
  iEvent.getByLabel(m_rctInput_tag,inputEmCands);
  
  Handle<L1CaloRegionCollection> inputRegions;
  iEvent.getByLabel(m_rctInput_tag,inputRegions);
  
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

  analyzeJets(cJets1, cJets2);
  analyzeJets(tJets1, tJets2);
  analyzeJets(fJets1, fJets2);

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
  float etMAng1 = etMPhi1*M_PI/36.;
  float etMAng2 = etMPhi2*M_PI/36.;

  float exMissRct = 0.0;
  float eyMissRct = 0.0;
  std::vector<float> exRctEtaBin(22,0.0);
  std::vector<float> eyRctEtaBin(22,0.0);
  for (L1CaloRegionCollection::const_iterator rgn = (*inputRegions).begin();
                                              rgn!= (*inputRegions).end(); rgn++) {
    float etRgn = static_cast<float>(rgn->et());
    float phBin = static_cast<float>(rgn->id().iphi());
    float phRgn = phBin*M_PI/9.;
    exMissRct -= etRgn*cos(phRgn);
    eyMissRct -= etRgn*sin(phRgn);

    unsigned etaBin = rgn->id().ieta();
    exRctEtaBin.at(etaBin) -= etRgn*cos(phRgn);
    eyRctEtaBin.at(etaBin) -= etRgn*sin(phRgn);
  }
  theMissEtVecRCT->Fill(exMissRct, eyMissRct);
  for (unsigned e=0; e<22; e++) {
    theRctExPerEtaBin.at(e)->Fill(exRctEtaBin.at(e));
    theRctEyPerEtaBin.at(e)->Fill(eyRctEtaBin.at(e));
  }

  theSumEtComparison->Fill ( etTot1,  etTot2 );
  theSumEtDifference->Fill ( etTot1 - etTot2 );
  theSumHtComparison->Fill ( etHad1,  etHad2 );
  theSumHtDifference->Fill ( etHad1 - etHad2 );
  theMissEtComparison->Fill(etMiss1,  etMiss2);
  theMissEtDifference->Fill(etMiss1 - etMiss2);
  theMissEtPhiGCT1->Fill(etMPhi1);
  theMissEtPhiGCT2->Fill(etMPhi2);
  theMissEtVecGCT1->Fill(etMiss1*cos(etMAng1), etMiss1*sin(etMAng1));
  theMissEtVecGCT2->Fill(etMiss2*cos(etMAng2), etMiss2*sin(etMAng2));
  theMETOffsetGCT1->Fill(etMiss1*cos(etMAng1)-exMissRct, etMiss1*sin(etMAng1)-eyMissRct);
  theMETOffsetGCT2->Fill(etMiss2*cos(etMAng2)-exMissRct, etMiss2*sin(etMAng2)-eyMissRct);

  theSumHtRatio->Fill(etHad1/etHad2);

}

void 
L1GctCompare::analyzeJets(const edm::Handle<L1GctJetCandCollection>& jColl1,
                          const edm::Handle<L1GctJetCandCollection>& jColl2) {

  L1GctJetCandCollection j1Matched;
  L1GctJetCandCollection j2Matched;
  L1GctJetCandCollection j1NonMatched;
  L1GctJetCandCollection j2NonMatched;

  L1GctJetCandCollection::const_iterator j1;
  L1GctJetCandCollection::const_iterator j2;

  static const uint16_t etaphiMask = 0x7fc0;

  j1=jColl1->begin();
  while (j1 != jColl1->end()) {
    if (j1->rank()==0) break;
    uint16_t pos1 = j1->raw() & etaphiMask;
    j2=jColl2->begin();
    bool match = false;
    while (j2 != jColl2->end()) {
      if (j2->rank()==0) break;
      uint16_t pos2 = j2->raw() & etaphiMask;
      if ( pos1==pos2 ) { match = true; break; }
      j2++;
    }
    if (match) j1Matched.push_back(*j1++);
    else j1NonMatched.push_back(*j1++);
  }

  j2=jColl2->begin();
  while (j2 != jColl2->end()) {
    if (j2->rank()==0) break;
    uint16_t pos2 = j2->raw() & etaphiMask;
    j1=jColl1->begin();
    bool match = false;
    while (j1 != jColl1->end()) {
      if (j1->rank()==0) break;
      uint16_t pos1 = j1->raw() & etaphiMask;
      if ( pos1==pos2 ) { match = true; break; }
      j1++;
    }
    if (match) j2Matched.push_back(*j2++);
    else j2NonMatched.push_back(*j2++);
  }
    
  // Now fill the histograms. First for matched pairs.
  assert (j1Matched.size()==j2Matched.size());
  j1 = j1Matched.begin();
  j2 = j2Matched.begin();
  while (j1 != j1Matched.end()) {
    int rcteta = j1->regionId().rctEta();
    int gcteta = (j1->etaSign()==0) ? (11 + rcteta) : (10 - rcteta);
    int gctphi = j1->phiIndex();
    // We have found a matched pair, fill some histograms
    float r1 = static_cast<float>(j1->rank());
    float r2 = static_cast<float>(j2->rank());
    theJetRankComparison->Fill(r1,r2);
    theJetRankDifference->Fill((r1-r2),r1);

    if (j1->isCentral()) {
      theCenJetRankComparison->Fill(r1,r2);
      theCenJetRankDifference->Fill((r1-r2),r1);
    }
    if (j1->isTau()) {
      theTauJetRankComparison->Fill(r1,r2);
      theTauJetRankDifference->Fill((r1-r2),r1);
    }
    if (j1->isForward()) {
      theFwdJetRankComparison->Fill(r1,r2);
      theFwdJetRankDifference->Fill((r1-r2),r1);
    }

    theJetRankDifferencePerEtaBin.at(rcteta)->Fill((r1-r2),r1);
    theEtaPhiMatchedJets->Fill(gcteta, gctphi);
    m_jets      += 1.0;
    m_totDiff   += (r1-r2);
    m_totDiffSq += (r1-r2)*(r1-r2);

    j1++;
    j2++;
  }
  j1 = j1NonMatched.begin();
  j2 = j2NonMatched.begin();
  while (j1 != j1NonMatched.end()) {
    int rcteta = j1->regionId().rctEta();
    int gcteta = (j1->etaSign()==0) ? (11 + rcteta) : (10 - rcteta);
    int gctphi = j1->phiIndex();
    theJetRankGCT1Only->Fill(j1->rank());
    theEtaPhiGCT1Only->Fill(gcteta, gctphi);
    j1++;
  }
  while (j2 != j2NonMatched.end()) {
    int rcteta = j2->regionId().rctEta();
    int gcteta = (j2->etaSign()==0) ? (11 + rcteta) : (10 - rcteta);
    int gctphi = j2->phiIndex();
    theJetRankGCT2Only->Fill(j2->rank());
    theEtaPhiGCT2Only->Fill(gcteta, gctphi);
    j2++;
  }
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

  theEtaPhiMatchedJets = dir.make<TH2F>("EtaPhiMatchedJets",  "Eta & phi for matched jets",
                                              22, 0., 22., 18, 0., 18.);
  theJetRankGCT1Only = dir.make<TH1F>("JetRankGCT1Only", "Jet Rank gct1, gct2 finds no jet",
                                              64, 0., 64.);
  theJetRankGCT2Only = dir.make<TH1F>("JetRankGCT2Only", "Jet Rank gct2, gct1 finds no jet",
                                              64, 0., 64.);
  theEtaPhiGCT1Only  = dir.make<TH2F>("EtaPhiGCT1Only",  "Eta & phi gct1, gct2 finds no jet",
                                              22, 0., 22., 18, 0., 18.);
  theEtaPhiGCT2Only  = dir.make<TH2F>("EtaPhiGCT2Only",  "Eta & phi gct2, gct1 finds no jet",
                                              22, 0., 22., 18, 0., 18.);

  

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
                                       128, 0.,512., 128, 0., 512.);
  theMissEtComparison = dir2.make<TH2F>("MissEtComparison", "Missing Et gct1 vs gct2",
                                       128, 0., 256., 128, 0., 256.);

  theSumEtDifference  = dir2.make<TH1F>("SumEtDifference",    "Total Et gct1-gct2", 100, -50., 50.);
  theSumHtDifference  = dir2.make<TH1F>("SumHtDifference",    "Total Ht gct1-gct2", 100, -50., 50.);
  theMissEtDifference = dir2.make<TH1F>("MissEtDifference", "Missing Et gct1-gct2", 100, -50., 50.);

  theMissEtPhiGCT1 = dir2.make<TH1F>("MissEtPhiGCT1", "Missing Et angle gct1", 72, 0., 72.);
  theMissEtPhiGCT2 = dir2.make<TH1F>("MissEtPhiGCT2", "Missing Et angle gct2", 72, 0., 72.);

  theMissEtVecRCT  = dir2.make<TH2F>("MissEtVecRCT", "Missing Et vector rct input regions",
                                     128, -256., 256., 128, -256., 256.);
  theMissEtVecGCT1 = dir2.make<TH2F>("MissEtVecGCT1", "Missing Et vector gct1",
                                     128, -256., 256., 128, -256., 256.);
  theMissEtVecGCT2 = dir2.make<TH2F>("MissEtVecGCT2", "Missing Et vector gct2",
                                     128, -256., 256., 128, -256., 256.);
  theMETOffsetGCT1 = dir2.make<TH2F>("METOffsetGCT1", "Missing Et gct1-rct",
                                     100, -20., 20., 100, -20., 20.);
  theMETOffsetGCT2 = dir2.make<TH2F>("METOffsetGCT2", "Missing Et gct2-rct",
                                     100, -20., 20., 100, -20., 20.);

  for (unsigned eta=0; eta<22; eta++) {
    std::stringstream ss;
    std::string title;
    ss << "RctExEtaBin" << eta;
    ss >> title;
    theRctExPerEtaBin.push_back(dir2.make<TH1F>(title.c_str(), "Rct Ex single eta bin", 128, -256., 256.));
  }
  for (unsigned eta=0; eta<22; eta++) {
    std::stringstream ss;
    std::string title;
    ss << "RctEyEtaBin" << eta;
    ss >> title;
    theRctEyPerEtaBin.push_back(dir2.make<TH1F>(title.c_str(), "Rct Ey single eta bin", 128, -256., 256.));
  }


  theSumHtRatio       = dir2.make<TH1F>("SumHtRatio",         "Total Ht gct1/gct2", 200, 0., 2.);
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
