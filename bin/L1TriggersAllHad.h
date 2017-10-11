#ifndef Analysis_TriggerStudies_L1TriggersAllHad_h
#define Analysis_TriggerStudies_L1TriggersAllHad_h 1

#include "Analysis/Core/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::core;


bool L1SingleJet(Analysis & , std::map<std::string, TH1F*> & , const int &, const float & );
bool L1DoubleJetC100(Analysis & , std::map<std::string, TH1F*> & , const std::string &, const float & );
bool L1DoubleJetC112(Analysis & , std::map<std::string, TH1F*> & , const std::string &, const float & );
bool L1DoubleJet100Eta2p4(Analysis & );
bool L1DoubleJet100Eta2p3_dEtaMax1p6(Analysis & );
bool L1DoubleJet112Eta2p3_dEtaMax1p6(Analysis & );
bool L1DoubleJet100Eta2p3_dEtaMax1p6_NotIn112(Analysis & );


// L1 triggers


bool L1SingleJet(Analysis & analysis, std::map<std::string, TH1F*> & hist, const int & ptmin , const float & pu)
{
   if ( ! analysis.triggerResult("HLT_L1SingleJet20_v") ) return false;
   std::string name  = "L1_SingleJet"+std::to_string(ptmin);
   
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
      if ( jet.pt() < float(ptmin) ) continue;
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 1 ) return false;
   
   if ( int(pu) == 30 ) 
   {
   
       hist[Form("jet1PT_%s",name.c_str())]  -> Fill(l1jet.at(0).pt());
       hist[Form("jet1ETA_%s",name.c_str())] -> Fill(l1jet.at(0).eta());
       hist[Form("jet1PHI_%s",name.c_str())] -> Fill(l1jet.at(0).phi());
      
      hist[Form("jetsN_%s",name.c_str())] -> Fill(l1jet.size());
      for ( size_t j = 0 ; j < l1jet.size() ; ++j )
      {
          TriggerObject jet = l1jet.at(j);
          hist[Form("jetsPT_%s",name.c_str())]  -> Fill(jet.pt());
          hist[Form("jetsETA_%s",name.c_str())] -> Fill(jet.eta());
          hist[Form("jetsPHI_%s",name.c_str())] -> Fill(jet.phi());
      }
   }
   
   return true;
   
}


// ----------------------------------------------------------------------
bool L1DoubleJetC100(Analysis & analysis, std::map<std::string, TH1F*> & hist, const std::string & name, const float & pu )
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   if ( int(pu) == 30 ) 
   {
       hist[Form("jet1PT_%s",name.c_str())]  -> Fill(l1jet.at(0).pt());
       hist[Form("jet1ETA_%s",name.c_str())] -> Fill(l1jet.at(0).eta());
       hist[Form("jet1PHI_%s",name.c_str())] -> Fill(l1jet.at(0).phi());
      
      hist[Form("jetsN_%s",name.c_str())] -> Fill(l1jet.size());
      for ( size_t j = 0 ; j < l1jet.size() ; ++j )
      {
          TriggerObject jet = l1jet.at(j);
          hist[Form("jetsPT_%s",name.c_str())]  -> Fill(jet.pt());
          hist[Form("jetsETA_%s",name.c_str())] -> Fill(jet.eta());
          hist[Form("jetsPHI_%s",name.c_str())] -> Fill(jet.phi());
      }
   }  
   
   return true;
}
bool L1DoubleJetC112(Analysis & analysis, std::map<std::string, TH1F*> & hist, const std::string & name, const float & pu )
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   std::vector<TriggerObject> l1jet2;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 112. ) continue;
      l1jet2.push_back(jet);
   }
   
   if ( l1jet2.size() < 2 )                          return false;
   
   
   if ( int(pu) == 30 ) 
   {
       hist[Form("jet1PT_%s",name.c_str())]  -> Fill(l1jet2.at(0).pt());
       hist[Form("jet1ETA_%s",name.c_str())] -> Fill(l1jet2.at(0).eta());
       hist[Form("jet1PHI_%s",name.c_str())] -> Fill(l1jet2.at(0).phi());
      
      hist[Form("jetsN_%s",name.c_str())] -> Fill(l1jet2.size());
      for ( size_t j = 0 ; j < l1jet2.size() ; ++j )
      {
          TriggerObject jet = l1jet2.at(j);
          hist[Form("jetsPT_%s",name.c_str())]  -> Fill(jet.pt());
          hist[Form("jetsETA_%s",name.c_str())] -> Fill(jet.eta());
          hist[Form("jetsPHI_%s",name.c_str())] -> Fill(jet.phi());
      }
   }  
   
   return true;
}

// ----------------------------------------------------------------------

bool L1DoubleJet100Eta2p4(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   std::vector<TriggerObject> l1jet2;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 100. || fabs(jet.eta()) > 2.4 ) continue;
      l1jet2.push_back(jet);
   }
   
   if ( l1jet2.size() < 2 )                          return false;
   
   
   return true;
   
}

// ----------------------------------------------------------------------

bool L1DoubleJet100Eta2p3_dEtaMax1p6(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   std::vector<TriggerObject> l1jet2;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 100. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet2.push_back(jet);
   }
   
   if ( l1jet2.size() < 2 )                          return false;
   
   bool deta = false;
   for ( int j1 = 0 ; j1 < (int)l1jet2.size()-1 ; ++j1 )
   {
      for ( int j2 = j1+1 ; j2 < (int)l1jet2.size() ; ++j2 )
      {
         float eta1 = l1jet2.at(j1).eta();
         float eta2 = l1jet2.at(j2).eta();
         if ( fabs(eta1-eta2) <= 1.6 ) deta = true;
      }
   }
   
   if ( ! deta )                                            return false;
   
   return true;
   
}
// ----------------------------------------------------------------------

bool L1DoubleJet112Eta2p3_dEtaMax1p6(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   std::vector<TriggerObject> l1jet2;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 112. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet2.push_back(jet);
   }
   
   if ( l1jet2.size() < 2 )                          return false;
   
   bool deta = false;
   for ( int j1 = 0 ; j1 < (int)l1jet2.size()-1 ; ++j1 )
   {
      for ( int j2 = j1+1 ; j2 < (int)l1jet2.size() ; ++j2 )
      {
         float eta1 = l1jet2.at(j1).eta();
         float eta2 = l1jet2.at(j2).eta();
         if ( fabs(eta1-eta2) <= 1.6 ) deta = true;
      }
   }
   
   if ( ! deta )                                            return false;
   
   return true;
   
}

bool L1DoubleJet100Eta2p3_dEtaMax1p6_NotIn112(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1DoubleJetC100_v") ) return false;
   
   auto l1dijet100 = analysis.collection<TriggerObject>("hltL1sDoubleJetC100");
   std::vector<TriggerObject> l1jet;
   
   TriggerObject prevJet;  // previous jet
   for ( int j = 0 ; j < l1dijet100->size() ; ++j )
   {
      TriggerObject jet = l1dijet100->at(j);
      if ( j > 0 )  // find possible duplicated object and skip it
      {
         prevJet = l1dijet100->at(j-1);
         if ( fabs(jet.pt()-prevJet.pt()) < 1.e-5 && fabs(jet.eta()-prevJet.eta()) < 1.e-5 && fabs(jet.phi()-prevJet.phi()) < 1.e-5 ) continue;
      }
      l1jet.push_back(jet);
   }
   
   if ( l1jet.size() < 2 ) return false;
   
   std::vector<TriggerObject> l1jet2;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 100. || fabs(jet.eta()) > 2.3 ) continue;
      l1jet2.push_back(jet);
   }
   
   if ( l1jet2.size() < 2 )                          return false;
   
   bool deta = false;
   for ( int j1 = 0 ; j1 < (int)l1jet2.size()-1 ; ++j1 )
   {
      for ( int j2 = j1+1 ; j2 < (int)l1jet2.size() ; ++j2 )
      {
         float eta1 = l1jet2.at(j1).eta();
         float eta2 = l1jet2.at(j2).eta();
         if ( fabs(eta1-eta2) <= 1.6 ) deta = true;
      }
   }
   
   if ( ! deta )                                            return false;
   
   // The new seed fired!
   // Check pure rate wrt L1_DoubleJetC112, if this also fired then do not count
   
   std::vector<TriggerObject> l1jet3;
   for ( size_t j = 0 ; j < l1jet.size() ; ++j )
   {
      TriggerObject jet = l1jet.at(j);
      if ( jet.pt() < 112. ) continue;
      l1jet3.push_back(jet);
   }
   
   if ( l1jet3.size() < 2 )                          return false;

   
   return true;
   
}


#endif  // Analysis_TriggerStudies_L1TriggersAllHad_h
