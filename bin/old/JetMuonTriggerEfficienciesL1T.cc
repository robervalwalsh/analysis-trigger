#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/Tools/bin/macro_config.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;

std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx);
std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx, const int * const nmin, const std::vector<float> * const ptmin, const std::vector<float> * const etamax );
bool MatchOnlineOffline(const Muon & muon, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Muon & muon, const std::vector<L1TMuon *> & tos);
bool MatchOnlineOffline(const Muon & muon, const std::shared_ptr< Collection<TriggerObject> > & tos);

bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::shared_ptr< Collection<TriggerObject> > & tos);

bool IsGoodPrimaryVertex(const Vertex & pv);

bool TriggerAccept(Analysis & analysis, const Jet & jet, const Muon & muon, const std::string & type = "nom");

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Primary vertex
   analysis.addTree<Vertex> ("PrimaryVertices",Form("%s/offlineSlimmedPrimaryVertices",treePath_.c_str()));
   // Jets
   analysis.addTree<Jet> ("Jets",Form("%s/%s",treePath_.c_str(),jetsCol_.c_str()));
   // L1T Jets
   analysis.addTree<L1TJet> ("L1TJets",Form("%s/%s",treePath_.c_str(),l1tjetsCol_.c_str()));
   // Muons
   analysis.addTree<Muon> ("Muons",Form("%s/%s",treePath_.c_str(),muonsCol_.c_str()));
   // L1T Muons
   analysis.addTree<L1TMuon> ("L1TMuons",Form("%s/%s",treePath_.c_str(),l1tmuonsCol_.c_str()));
   // Trigger path info
   analysis.triggerResults(Form("%s/%s",treePath_.c_str(),triggerCol_.c_str()));
   // Trigger objects
   for ( auto & obj : triggerObjectsJets_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   }
   for ( auto & obj : triggerObjectsMuons_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   }
   
   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   
   h1["dR_num"]      = new TH1F("dR_num"    , "" , 50,  0   ,  1. );
   h1["dR_den"]      = new TH1F("dR_den"    , "" , 50,  0   ,  1. );
   h1["deta_num"]    = new TH1F("deta_num"  , "" , 50,  0   ,  1. );
   h1["deta_den"]    = new TH1F("deta_den"  , "" , 50,  0   ,  1. );
   h1["dphi_num"]    = new TH1F("dphi_num"  , "" , 50,  0   ,  1. );
   h1["dphi_den"]    = new TH1F("dphi_den"  , "" , 50,  0   ,  1. );
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%10000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      
// Offline selection
// =================
      // Object - std::shared_ptr< Collection<Object> >
      auto pvs = analysis.collection<Vertex>("PrimaryVertices");
      if ( pvs->size() < 1 ) continue;
      if ( ! IsGoodPrimaryVertex(pvs->at(0)) ) continue;
      
      // Muons
      auto slimmedMuons = analysis.collection<Muon>("Muons");      
      // selection of muons
      std::vector<Muon *> selectedIdMuons;
      for ( int m = 0 ; m < slimmedMuons->size() ; ++m )
      {
         Muon * muon = &(slimmedMuons->at(m));
         if ( muon -> isMediumMuon() )
            selectedIdMuons.push_back(muon);
      }
      if ( (int)selectedIdMuons.size() < nmuonsmin_ ) continue;
      
      std::vector<Muon *> selectedMuons;
      // selection only on the leading  id'ed muon
      for ( int m = 0 ; m < (int)selectedIdMuons.size() ; ++m )
      {
         Muon * muon = selectedIdMuons.at(m);
         if ( muon->pt() > muonsptmin_[0] && muon->pt() <= muonsptmax_[0] && fabs(muon->eta()) < muonsetamax_[0] )
         {
            selectedMuons.push_back(muon);
         }
      }
      if ( (int)selectedMuons.size() < nmuonsmin_ ) continue;
      
      // Jets
      auto slimmedJets = analysis.collection<Jet>("Jets");      
      // selection of jets
      std::vector<Jet *> selectedLooseJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         Jet * jet = &(slimmedJets->at(j));
         if ( jet -> idLoose() )
            selectedLooseJets.push_back(jet);
      }
      if ( (int)selectedLooseJets.size() < njetsmin_ ) continue;
      std::vector<Jet *> selectedJets;
      // selection only on the leading loose id'ed jets
      for ( int j = 0 ; j < (int)selectedLooseJets.size() ; ++j )
      {
         Jet * jet = selectedLooseJets.at(j);
         if ( jet->pt() > jetsptmin_[0] && jet->pt() <= jetsptmax_[0] && fabs(jet->eta()) < jetsetamax_[0] )
            selectedJets.push_back(jet);
      }
      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
      Jet * jet = selectedJets.at(0);
      Muon * muon = selectedMuons.at(0);
      
      float deltar = -1;
      for ( auto & m : selectedMuons )
      {
          deltar = m->deltaR(*jet);
          if ( deltar < drmax_ )
          {
             muon = m;
             break;
          }
      }
      if ( deltar < 0 || deltar > drmax_ ) continue;
      
// Trigger analysis, reference trigger (denominator)
// =================================================  
      bool refAccept = TriggerAccept(analysis,*jet, *muon,"ref");
      if ( ! refAccept ) continue;
      h1["dR_den"] -> Fill(deltar);
      
      bool nomAccept = TriggerAccept(analysis,*jet, *muon,"nom");
      if ( ! nomAccept ) continue;
      h1["dR_num"] -> Fill(deltar);
      
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
   std::cout << "JetMuonTriggerEfficiencies: program finished" << std::endl;
   
} //end main


// ================================

bool TriggerAccept(Analysis & analysis, const Jet & jet, const Muon & muon, const std::string & type)
{
   
// Trigger analysis, reference trigger (denominator)
// =================================================      
   bool hltAccept = analysis.triggerResult(hltPaths_[0]) && analysis.triggerResult(hltPaths_[1]);
   if ( !hltAccept ) return false; 
   
   // L1 objects
   // Emulation with L1T candidate --- L1TMuons
   std::vector<L1TMuon *> l1muons;
   if ( l1tmuonsnmin_ > 0 )
   {
      auto l1tmuons = analysis.collection<L1TMuon>("L1TMuons");
      for ( int j = 0 ; j < l1tmuons->size() ; ++j )
      {
         L1TMuon * l1tmuon = &(l1tmuons->at(j));
         if ( l1tmuon->pt() >= l1tmuonsptmin_[0] && fabs(l1tmuon->eta()) <= l1tmuonsetamax_[0] )
            l1muons.push_back(l1tmuon);
      }
      if ( (int)l1muons.size() < l1tmuonsnmin_ ) return false;
   }
   
   // Emulation with L1T candidate --- L1TJets
   std::vector<L1TJet *> l1jets;
   if ( l1tjetsnmin_ > 0 )
   {
      auto l1tjets = analysis.collection<L1TJet>("L1TJets");
      for ( int j = 0 ; j < l1tjets->size() ; ++j )
      {
         L1TJet * l1tjet = &(l1tjets->at(j));
         if ( l1tjet->pt() >= l1tjetsptmin_[0] && fabs(l1tjet->eta()) <= l1tjetsetamax_[0] )
            l1jets.push_back(l1tjet);
      }
      if ( (int)l1jets.size() < l1tjetsnmin_ ) return false;
   }
   
   // L3 objects
   auto l3muons = analysis.collection<TriggerObject>(triggerObjectsMuons_[0]);
   auto l3jets = analysis.collection<TriggerObject>(triggerObjectsJets_[0]);
   
   // match leading muon to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
   if ( matchonoff_ )
   {
      // matching muons
      if (  !(MatchOnlineOffline(muon,l1muons) && MatchOnlineOffline(muon,l3muons)) ) return false;
      // matching jets
      if (  !(MatchOnlineOffline(jet,l1jets) && MatchOnlineOffline(jet,l3jets)) ) return false;
      
   }
   
   if ( type == "ref" ) return true;
   
   // emulate delta R between jet-muon at trigger level - at least ONE muon-jet
   // L1
   bool accept = false;
   for ( auto m : l1muons )
   {
      for ( auto j : l1jets )
      {
         if ( m->deltaR(*j) <= 0.4 ) accept = true ;
      }
   }
   
   if ( ! accept ) return false;
   
   // L3
   accept = false;
   for ( int m = 0 ; m < l3muons->size() ; ++m )
   {
      TriggerObject mu = l3muons->at(m);
      for ( int j = 0 ; j < l3jets->size() ; ++j )
      {
         TriggerObject je = l3jets->at(j);
         if ( mu.deltaR(je) <= 0.5 ) accept = true ; // 0.5 is the threshold of the BTagMu module
      }
   }
   
   if ( ! accept ) return false;
   
   return true; 

}



std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx)
{
   std::vector<TriggerObject *> seltos;
   
   // select trigger objects
   for ( int j = 0; j < tos->size() ; ++j )
   {
      TriggerObject * to = &tos->at(j);
      if ( tonmin_[idx] < 1 )
      {
         seltos.push_back(to);
      }
      else  // emulation
      {
         if ( to -> pt() >= toptmin_[idx][0] && fabs( to -> eta() ) <= toetamax_[idx][0] ) 
            seltos.push_back(to);
      }
   }
   return seltos;

}

std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx, const int * const nmin, const std::vector<float> * const ptmin, const std::vector<float> * const etamax )
{
   std::vector<TriggerObject *> seltos;
   
   // select trigger objects
   for ( int j = 0; j < tos->size() ; ++j )
   {
      TriggerObject * to = &tos->at(j);
      if ( nmin[idx] < 1 )
      {
         seltos.push_back(to);
      }
      else  // emulation
      {
         if ( to -> pt() >= ptmin[idx][0] && fabs( to -> eta() ) <= etamax[idx][0] ) 
            seltos.push_back(to);
      }
   }
   return seltos;

}

bool MatchOnlineOffline(const Muon & muon, const std::vector<TriggerObject *> & tos)
{
   bool match = false;
   if ( tos.size() < 1 ) return true; // ?????
   
   for ( size_t i = 0; i < tos.size() ; ++i )
   {
      if ( muon.deltaR(*tos[i]) < matchonoffdrmax_ ) match = true;
      if ( match ) break;
   }
   return match;
   
}

bool MatchOnlineOffline(const Muon & cand, const std::shared_ptr< Collection<TriggerObject> > & tos)
{
   bool match = false;
   if ( tos->size() < 1 ) return true; // ?????
   
   for ( int i = 0; i < tos->size() ; ++i )
   {
      if ( cand.deltaR(tos->at(i)) < matchonoffdrmax_ ) match = true;
      if ( match ) break;
   }
   return match;
   
}

bool MatchOnlineOffline(const Jet & cand, const std::shared_ptr< Collection<TriggerObject> > & tos)
{
   bool match = false;
   if ( tos->size() < 1 ) return true; // ?????
   
   for ( int i = 0; i < tos->size() ; ++i )
   {
      if ( cand.deltaR(tos->at(i)) < matchonoffdrmax_ ) match = true;
      if ( match ) break;
   }
   return match;
   
}

bool MatchOnlineOffline(const Muon & muon, const std::vector<L1TMuon *> & tos)
{
   bool match = false;
   if ( tos.size() < 1 ) return true; // ?????
   
   for ( size_t i = 0; i < tos.size() ; ++i )
   {
      if ( muon.deltaR(*tos[i]) < matchonoffdrmax_ ) match = true;
      if ( match ) break;
   }
   return match;
   
}

bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos)
{
   bool match = false;
   if ( tos.size() < 1 ) return true; // ?????
   
   for ( size_t i = 0; i < tos.size() ; ++i )
   {
      if ( jet.deltaR(*tos[i]) < matchonoffdrmax_ ) match = true;
      if ( match ) break;
   }
   return match;
   
}

bool IsGoodPrimaryVertex(const Vertex & pv)
{
   bool good = false;
   
   if ( ! pv.fake() && pv.ndof() > 4 && fabs(pv.z()) <= 24 && pv.rho() <=2 ) good = true;
   else good = false;
   
   return good;
}
