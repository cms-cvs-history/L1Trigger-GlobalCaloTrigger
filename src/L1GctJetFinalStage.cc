#include "L1Trigger/GlobalCaloTrigger/interface/L1GctJetFinalStage.h"

L1GctJetFinalStage::L1GctJetFinalStage()
{
}

L1GctJetFinalStage::~L1GctJetFinalStage()
{
}

void L1GctJetFinalStage::reset()
{
}

void L1GctJetFinalStage::process()
{
}

void L1GctJetFinalStage::setInputJet(int i, L1GctJet jet)
{
    inputJets.push_back(jet);
}
