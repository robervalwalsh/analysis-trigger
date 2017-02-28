#ifndef Analysis_TriggerStudies_L1TriggersAllHad_h
#define Analysis_TriggerStudies_L1TriggersAllHad_h 1

#include "Analysis/Core/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


bool L1DoubleJetC100(Analysis & , std::map<std::string, TH1F*> & , const std::string & );
bool L1DoubleJet100Eta2p3(Analysis & );
bool L1DoubleJet100Eta2p3_dEtaMax1p6(Analysis & );


// L1 triggers

// ----------------------------------------------------------------------
bool L1DoubleJetC100(Analysis & analysis, std::map<std::string, TH1F*> & hist, const std::string & name )
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   
   
   hist[Form("jetsN_%s",name.c_str())] -> Fill(l1dijet100->size());
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
       TriggerObject jet = l1dijet100->at(j);
       hist[Form("jetsPT_%s",name.c_str())]  -> Fill(jet.pt());
       hist[Form("jetsETA_%s",name.c_str())] -> Fill(jet.eta());
       hist[Form("jetsPHI_%s",name.c_str())] -> Fill(jet.phi());
   }
   
   
   return true;
}

// ----------------------------------------------------------------------

bool L1DoubleJet100Eta2p3(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   if ( l1dijet100->size() < 2 )                           return false;
   
   std::vector<TriggerObject> l1dijetc100v2;
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( jet.pt() < 100. || fabs(jet.eta()) > 2.3 ) continue;
      l1dijetc100v2.push_back(jet);
   }
   
   if ( l1dijetc100v2.size() < 2 )                          return false;
   
   return true;
   
}

// ----------------------------------------------------------------------

bool L1DoubleJet100Eta2p3_dEtaMax1p6(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   if ( l1dijet100->size() < 2 )                            return false;
   
   std::vector<TriggerObject> l1dijetc100v2;
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( jet.pt() < 100. || fabs(jet.eta()) > 2.3 ) continue;
      l1dijetc100v2.push_back(jet);
   }
   
   if ( l1dijetc100v2.size() < 2 )                          return false;
   
   bool deta = false;
   for ( int j1 = 0 ; j1 < (int)l1dijetc100v2.size()-1 ; ++j1 )
   {
      for ( int j2 = j1+1 ; j2 < (int)l1dijetc100v2.size() ; ++j2 )
      {
         float eta1 = l1dijetc100v2.at(j1).eta();
         float eta2 = l1dijetc100v2.at(j2).eta();
         if ( fabs(eta1-eta2) <= 1.6 ) deta = true;
      }
   }
   
   if ( ! deta )                                            return false;
   
   return true;
   
}
#endif  // Analysis_TriggerStudies_L1TriggersAllHad_h
