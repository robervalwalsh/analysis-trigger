#ifndef Analysis_TriggerStudies_L1TriggersSemiLep_h
#define Analysis_TriggerStudies_L1TriggersSemiLep_h 1

#include "Analysis/Tools/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


bool L1SingleMu3(Analysis & );
bool L1Mu3JetC60(Analysis & );
bool L1Mu10DiJet32(Analysis & );
bool L1Mu12DiJet32(Analysis & );
bool L1Mu10DiJet40(Analysis & );


// L1 triggers

// ----------------------------------------------------------------------
bool L1SingleMu3(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1SingleMu3_v") ) return false;
   
   return true;
}

bool L1Mu10DiJet32(Analysis & analysis )
{
   if ( ! ( analysis.triggerResult("HLT_L1SingleMu3_v") || analysis.triggerResult("HLT_L1SingleJet20_v") ) ) return false;

   auto l1jet20 = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   
   std::vector<TriggerObject> l1jet;
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1jet20->size() ; ++j )
   {
      TriggerObject jet = l1jet20->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1jet20->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && jet.p4().DeltaPhi(prevJet.p4()) < 1.e-5 ) continue;
      }
      
      // kinematic selection
      if ( jet.pt() < 32. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   // muon in jet
   auto l1mu3 = analysis.collection<TriggerObject>("hltL1sSingleMu3");
   std::vector<TriggerObject> l1mu;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      for ( int m = 0; m < l1mu3->size() ; ++m )
      {
         TriggerObject mu = l1mu3->at(m);
         if ( mu.pt() < 10. || fabs(mu.eta()) > 2.3 ) continue;
         if ( jet.p4().DeltaR(mu.p4()) > 0.4 ) continue;
//         if ( fabs(jet.eta()-mu.eta()) > 0.4 || jet.p4().DeltaPhi(mu.p4()) > 0.4 ) continue;
         l1mu.push_back(mu);
      }
   }
   
   if ( l1mu.size() < 1 ) return false;
        
   return true;
}

bool L1Mu12DiJet32(Analysis & analysis )
{
   if ( ! ( analysis.triggerResult("HLT_L1SingleMu3_v") || analysis.triggerResult("HLT_L1SingleJet20_v") ) ) return false;

   auto l1jet20 = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   
   std::vector<TriggerObject> l1jet;
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1jet20->size() ; ++j )
   {
      TriggerObject jet = l1jet20->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1jet20->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && jet.p4().DeltaPhi(prevJet.p4()) < 1.e-5 ) continue;
      }
      
      // kinematic selection
      if ( jet.pt() < 32. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   // muon in jet
   auto l1mu3 = analysis.collection<TriggerObject>("hltL1sSingleMu3");
   std::vector<TriggerObject> l1mu;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      for ( int m = 0; m < l1mu3->size() ; ++m )
      {
         TriggerObject mu = l1mu3->at(m);
         if ( mu.pt() < 12. || fabs(mu.eta()) > 2.3 ) continue;
         if ( jet.p4().DeltaR(mu.p4()) > 0.4 ) continue;
//         if ( fabs(jet.eta()-mu.eta()) > 0.4 || jet.p4().DeltaPhi(mu.p4()) > 0.4 ) continue;
         l1mu.push_back(mu);
      }
   }
   
   if ( l1mu.size() < 1 ) return false;
        
   return true;
}

bool L1Mu10DiJet40(Analysis & analysis )
{
   if ( ! ( analysis.triggerResult("HLT_L1SingleMu3_v") || analysis.triggerResult("HLT_L1SingleJet20_v") ) ) return false;

   auto l1jet20 = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   
   std::vector<TriggerObject> l1jet;
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1jet20->size() ; ++j )
   {
      TriggerObject jet = l1jet20->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1jet20->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && jet.p4().DeltaPhi(prevJet.p4()) < 1.e-5 ) continue;
      }
      
      // kinematic selection
      if ( jet.pt() < 40. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   // muon in jet
   auto l1mu3 = analysis.collection<TriggerObject>("hltL1sSingleMu3");
   std::vector<TriggerObject> l1mu;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      for ( int m = 0; m < l1mu3->size() ; ++m )
      {
         TriggerObject mu = l1mu3->at(m);
         if ( mu.pt() < 10. || fabs(mu.eta()) > 2.3 ) continue;
         if ( jet.p4().DeltaR(mu.p4()) > 0.4 ) continue;
//         if ( fabs(jet.eta()-mu.eta()) > 0.4 || jet.p4().DeltaPhi(mu.p4()) > 0.4 ) continue;
         l1mu.push_back(mu);
      }
   }
   
   if ( l1mu.size() < 1 ) return false;
        
   return true;
}


bool L1Mu3JetC60(Analysis & analysis )
{
   if ( ! ( analysis.triggerResult("HLT_L1SingleMu3_v") || analysis.triggerResult("HLT_L1SingleJet20_v") ) ) return false;

   auto l1jet20 = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   
   std::vector<TriggerObject> l1jet;
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1jet20->size() ; ++j )
   {
      TriggerObject jet = l1jet20->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1jet20->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      
      // kinematic selection
      if ( jet.pt() < 60 || fabs(jet.eta()) > 3. ) continue;
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 1 ) return false;
   
   // muon in jet
   auto l1mu3 = analysis.collection<TriggerObject>("hltL1sSingleMu3");
   std::vector<TriggerObject> l1mu;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      for ( int m = 0; m < l1mu3->size() ; ++m )
      {
         TriggerObject mu = l1mu3->at(m);
         if ( fabs(jet.eta()-mu.eta()) > 0.4 || jet.p4().DeltaPhi(mu.p4()) > 0.4 ) continue;
         l1mu.push_back(mu);
      }
   }
   
   if ( l1mu.size() < 1 ) return false;
        
   return true;
}

#endif  // Analysis_TriggerStudies_L1TriggersSemiLep_h
