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
bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos);

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Jets
   analysis.addTree<Jet> ("Jets",Form("%s/%s",treePath_.c_str(),jetsCol_.c_str()));
   // L1T Jets
   analysis.addTree<L1TJet> ("L1TJets",Form("%s/%s",treePath_.c_str(),l1tjetsCol_.c_str()));
   // Trigger path info
   analysis.triggerResults(Form("%s/%s",treePath_.c_str(),triggerCol_.c_str()));
   // Trigger objects
   for ( auto & obj : triggerObjects_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   }
   
   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   
   h1["n"]   = new TH1F("n"   , "" ,   5,  0   ,   5  );
   h1["pt"]  = new TH1F("pt"  , "" , 500,  0   ,  1000  );
   h1["eta"] = new TH1F("eta" , "" , 100, -2.5 ,   2.5);
   h1["phi"] = new TH1F("phi" , "" , 128, -3.2 ,   3.2);
   
   if ( psweight_ )
   {
       h1["n_psw"]   = new TH1F("n_psw"   , "" ,   5,  0   ,   5  );
       h1["pt_psw"]  = new TH1F("pt_psw"  , "" , 500,  0   ,  1000  );
       h1["eta_psw"] = new TH1F("eta_psw" , "" , 100, -2.5 ,   2.5);
       h1["phi_psw"] = new TH1F("phi_psw" , "" , 128, -3.2 ,   3.2);
   
   }

      
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      
      if (! isMC_ )
      {
         int json_status = analysis.processJsonFile(json_);
         if ( json_status < 0 ) 
         {
            std::cout << "Error from processing json. Please check your json file." << std::endl;
            return -1;
         }
      }

// Offline selection
// =================
      // Object - std::shared_ptr< Collection<Object> >
      
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
      
      
// Trigger analysis
// ==================================      
      bool hltAccept = analysis.triggerResult(hltPath_);
      int hltPs = 1;
      
      if ( psweight_ ) hltPs = analysis.triggerPrescale(hltPath_);

      bool l1Accept  = true;
      int l1Ps = 1;
      if ( l1Seed_ != "" )
      {
         analysis.triggerResult(l1Seed_);
         if ( psweight_ ) l1Ps  = analysis.triggerPrescale(l1Seed_);
      }
      
      int psw = l1Ps * hltPs;
      
      if ( !hltAccept || !l1Accept ) continue; // Selecting one of the OR L1 seeds.
      
      
      // Emulation with L1T candidate
      // L1TJets
      // selection of jets
      std::vector<L1TJet *> selectedL1TJets;
      if ( l1tjetsnmin_ > 0 )
      {
         auto l1tjets = analysis.collection<L1TJet>("L1TJets");
         for ( int j = 0 ; j < l1tjets->size() ; ++j )
         {
            L1TJet * l1tjet = &(l1tjets->at(j));
            if ( l1tjet->pt() > l1tjetsptmin_[0] && fabs(l1tjet->eta()) < l1tjetsetamax_[0] )
               selectedL1TJets.push_back(l1tjet);
         }
         if ( (int)selectedL1TJets.size() < l1tjetsnmin_ ) continue;
      }
      
      // Emulation with trigger objects
      std::vector<TriggerObject *> selectedL1Jets;
      std::vector<TriggerObject *> selectedL2Jets;
      std::vector<TriggerObject *> selectedL3Jets;
      
      if ( triggerObjects_.size() > 2 )
      {
         // select L1 objects
         auto l1jets = analysis.collection<TriggerObject>(triggerObjects_[0]);
         selectedL1Jets = SelectTriggerObjects(l1jets,0);
         if ( (int)selectedL1Jets.size() < tonmin_[0]  ) continue;
         
         // select L2 objects
         auto l2jets = analysis.collection<TriggerObject>(triggerObjects_[1]);
         selectedL2Jets = SelectTriggerObjects(l2jets,1);
         if ( (int)selectedL2Jets.size() < tonmin_[1]  ) continue;
      
         // select L3 objects
         auto l3jets = analysis.collection<TriggerObject>(triggerObjects_[2]);
         selectedL3Jets = SelectTriggerObjects(l3jets,0);
         if ( (int)selectedL3Jets.size() < tonmin_[2]  ) continue;
      }
      
      // match leading jet to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
      if ( matchonoff_ )
      {
         bool l1match = false;
         if ( tonmin_[0] > 0 )       l1match = MatchOnlineOffline(*(jet),selectedL1Jets);
         else if (l1tjetsnmin_ > 0)  l1match = MatchOnlineOffline(*(jet),selectedL1TJets);
         else l1match = true;
         if ( !l1match ) continue;
         
         bool l2match = MatchOnlineOffline(*(jet),selectedL2Jets);
         if ( !l2match ) continue;
         
         bool l3match = MatchOnlineOffline(*(jet),selectedL3Jets);
         if ( !l3match ) continue;
      }
      
      // fill histograms
      h1["pt" ]     -> Fill(jet->pt() );
      h1["eta"]     -> Fill(jet->eta());
      h1["phi"]     -> Fill(jet->phi());
      if ( psweight_ )
      {
         h1["pt_psw" ] -> Fill(jet->pt() , psw);
         h1["eta_psw"] -> Fill(jet->eta(), psw);
         h1["phi_psw"] -> Fill(jet->phi(), psw);
      }
      
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
   std::cout << "SingleJetTriggerAnalysis: program finished" << std::endl;
   
} //end main


// ================================

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

bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos)
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
