#include "L1Trigger/GlobalCaloTrigger/plugins/L1GctValidation.h"

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/L1TObjects/interface/L1CaloEtScale.h"
#include "CondFormats/L1TObjects/interface/L1GctJetEtCalibrationFunction.h"
#include "CondFormats/DataRecord/interface/L1JetEtScaleRcd.h"
#include "CondFormats/DataRecord/interface/L1GctJetCalibFunRcd.h"

#include <math.h>

L1GctValidation::L1GctValidation(const edm::ParameterSet& iConfig) :
   m_energy_tag(iConfig.getUntrackedParameter<std::string>("L1GctEmulDigis","l1GctEmulDigis"))
{
}


L1GctValidation::~L1GctValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
L1GctValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   std::cout << "Hello from L1GctValidation::analyze()" << std::endl;

  // Get the scales from the event setup
  ESHandle< L1GctJetEtCalibrationFunction > calibFun ;
  iSetup.get< L1GctJetCalibFunRcd >().get( calibFun ) ; // which record?
  ESHandle< L1CaloEtScale > etScale ;
  iSetup.get< L1JetEtScaleRcd >().get( etScale ) ; // which record?

  double lsbForEt = etScale.product()->linearLsb();
  double lsbForHt = calibFun.product()->getHtScaleLSB();

  // Get the Gct energy sums from the event
  Handle< L1GctEtTotal > sumEt ;
  iEvent.getByLabel( m_energy_tag, sumEt ) ;
  Handle< L1GctEtHad > sumHt ;
  iEvent.getByLabel( m_energy_tag, sumHt ) ;
  Handle< L1GctEtMiss > missEt ;
  iEvent.getByLabel( m_energy_tag, missEt ) ;

  double etTot  = static_cast<double>(sumEt->et());
  double etHad  = static_cast<double>(sumHt->et());
  double etMiss = static_cast<double>(missEt->et());
  int phibin = missEt->phi();
  if (phibin>=36) phibin -= 72;
  double etMPhi = static_cast<double>(phibin);
  double etMAng = (etMPhi+0.5)*M_PI/36.;

  theSumEtInLsb->Fill(etTot);
  theSumHtInLsb->Fill(etHad);
  theMissEtInLsb->Fill(etMiss);
  theSumEtInGeV->Fill(etTot*lsbForEt);
  theSumHtInGeV->Fill(etHad*lsbForHt);
  theMissEtInGeV->Fill(etMiss*lsbForEt);
  theMissEtAngle->Fill(etMAng);
  theMissEtVector->Fill(etMiss*lsbForEt*cos(etMAng),etMiss*lsbForEt*sin(etMAng));

  // Get jet counts from the event
  Handle< L1GctJetCounts > jetCountDigi ;
  iEvent.getByLabel( m_energy_tag, jetCountDigi ) ;

  for (unsigned jc=0; jc<L1GctJetCounts::MAX_TOTAL_COUNTS; jc++) {
    theJetCounts.at(jc)->Fill(jetCountDigi->count(jc));
  }

  theHfEtSumPositiveEta->Fill(jetCountDigi->hfTowerEtSumPositiveEta());  
  theHfEtSumNegativeEta->Fill(jetCountDigi->hfTowerEtSumNegativeEta());  
  theHfTowerCountPositiveEta->Fill(jetCountDigi->hfTowerCountPositiveEta());  
  theHfTowerCountNegativeEta->Fill(jetCountDigi->hfTowerCountNegativeEta());  


}

// ------------ method called once each job just before starting event loop  ------------
void 
L1GctValidation::beginJob(const edm::EventSetup&)
{
  edm::Service<TFileService> fs;

  TFileDirectory dir0 = fs->mkdir("L1GctEtSums");

  theSumEtInLsb   = dir0.make<TH1F>("SumEtInLsb",   "Total Et (GCT units)",
                                    128, 0., 2048.); 
  theSumHtInLsb   = dir0.make<TH1F>("SumHtInLsb",   "Total Ht (GCT units)",
                                    128, 0., 2048.); 
  theMissEtInLsb  = dir0.make<TH1F>("MissEtInLsb",  "Missing Et magnitude (GCT units)",
                                    128, 0., 1024.); 
  theSumEtInGeV   = dir0.make<TH1F>("SumEtInGeV",   "Total Et (in GeV)",
                                    100, 0., 1000.); 
  theSumHtInGeV   = dir0.make<TH1F>("SumHtInGeV",   "Total Ht (in GeV)",
                                    100, 0., 1000.); 
  theMissEtInGeV  = dir0.make<TH1F>("MissEtInGeV",  "Missing Et magnitude (in GeV)",
                                    100, 0., 500.); 
  theMissEtAngle  = dir0.make<TH1F>("MissEtAngle",  "Missing Et angle",
                                    72, -M_PI, M_PI);
  theMissEtVector = dir0.make<TH2F>("MissEtVector", "Missing Ex vs Missing Ey",
                                    100, -100., 100., 100, -100., 100.); 

  TFileDirectory dir1 = fs->mkdir("L1GctJetCounts");

  for (unsigned jc=0; jc<L1GctJetCounts::MAX_TOTAL_COUNTS; jc++) {
    std::stringstream ss;
    std::string title;
    std::string header;
    ss << "JetCount#" << jc;
    ss >> title;
    ss << "Jet Count number " << jc;
    if (jc== 6 || jc== 7) { ss << " (Hf tower count)"; }
    if (jc==10 || jc==11) { ss << " (Hf Et sum MSB)"; }
    ss >> header;
    theJetCounts.push_back(dir1.make<TH1F>(title.c_str(), header.c_str(), 32, 0., 32.));
  }

  theHfEtSumPositiveEta      = dir1.make<TH1F>("HfEtSumPositiveEta", "Hf Inner Ring Et eta+",
                                               100, 0., 200.);
  theHfEtSumNegativeEta      = dir1.make<TH1F>("HfEtSumNegativeEta", "Hf Inner Ring Et eta-",
                                               100, 0., 200.);
  theHfTowerCountPositiveEta = dir1.make<TH1F>("HfTowerCountPositiveEta", "Hf Threshold bits eta+",
                                               20, 0., 20.);
  theHfTowerCountNegativeEta = dir1.make<TH1F>("HfTowerCountNegativeEta", "Hf Threshold bits eta-",
                                               20, 0., 20.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1GctValidation::endJob() {
}

DEFINE_ANOTHER_FWK_MODULE(L1GctValidation);

