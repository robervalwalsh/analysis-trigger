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
   h1["pt"]  = new TH1F("pt"  , "" , 100,  0   ,  1000  );
   h1["eta"] = new TH1F("eta" , "" ,  12, -2.4 ,   2.4);
   h1["phi"] = new TH1F("phi" , "" ,  16, -3.2 ,   3.2);
      
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
         if ( jet->pt() > jetsptmin_[0] && fabs(jet->eta()) < jetsetamax_[0] )
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
      auto l1tjets = analysis.collection<L1TJet>("L1TJets");
      // selection of jets
      std::vector<L1TJet *> selectedL1TJets;
      for ( int j = 0 ; j < l1tjets->size() ; ++j )
      {
         L1TJet * l1tjet = &(l1tjets->at(j));
         if ( l1tjet->pt() > l1tjetsptmin_[0] && fabs(l1tjet->eta()) < l1tjetsetamax_[0] )
            selectedL1TJets.push_back(l1tjet);
      }
      if ( (int)selectedL1TJets.size() < l1tjetsnmin_ ) continue;
      
      // Emulation with trigger objects
      std::vector<TriggerObject *> selectedL1Jets;
      std::vector<TriggerObject *> selectedL2Jets;
      std::vector<TriggerObject *> selectedL3Jets;
      
      if ( triggerObjects_.size() > 2 )
      {
         // select L1 objects
//         selectedL1Jets = SelectTriggerObjects(analysis,0);
      
         auto l1jets = analysis.collection<TriggerObject>(triggerObjects_[0]);
         for ( int j = 0; j < l1jets->size() ; ++j )
         {
            TriggerObject * l1jet = &l1jets->at(j);
            if ( ntomin_[0] < 1 )
            {
               selectedL1Jets.push_back(l1jet);
            }
            else  // L1 emulation
            {
               if ( l1jet -> pt() >= toptmin_[0][0] && fabs( l1jet -> eta() ) <= toetamax_[0][0] ) 
                  selectedL1Jets.push_back(l1jet);
            }
         }
         if ( (int)selectedL1Jets.size() < ntomin_[0]  ) continue;
         
         // select L2 objects
         auto l2jets = analysis.collection<TriggerObject>(triggerObjects_[1]);
         for ( int j = 0; j < l2jets->size() ; ++j )
         {
            TriggerObject * l2jet = &l2jets->at(j);
            if ( ntomin_[1] < 1 )
            {
               selectedL2Jets.push_back(l2jet);
            }
            else  // L2 emulation
            {
               if ( l2jet -> pt() >= toptmin_[1][0] && fabs( l2jet -> eta() ) <= toetamax_[1][0] ) 
                  selectedL2Jets.push_back(l2jet);
            }
         }
         if ( (int)selectedL2Jets.size() < ntomin_[1]  ) continue;
      
         // select L3 objects
         auto l3jets = analysis.collection<TriggerObject>(triggerObjects_[2]);
         for ( int j = 0; j < l3jets->size() ; ++j )
         {
            TriggerObject * l3jet = &l3jets->at(j);
            if ( ntomin_[2] < 1 )
            {
               selectedL3Jets.push_back(l3jet);
            }
            else  // L3 emulation
            {
               if ( l3jet -> pt() >= toptmin_[2][0] && fabs( l3jet -> eta() ) <= toetamax_[2][0] ) 
                  selectedL3Jets.push_back(l3jet);
            }
         }
         if ( (int)selectedL3Jets.size() < ntomin_[2]  ) continue;
      }
      
      // match leading jet to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
      if ( triggerObjects_.size() > 2 )
      {
         bool l1match = false;
         if ( l1tjetsnmin_ > 0 )
         {
            for ( size_t l1 = 0; l1 < selectedL1TJets.size() ; ++l1 )
            {
               if ( jet -> deltaR(*selectedL1TJets[l1]) < matchonoffdrmax_ ) l1match = true;
               if ( l1match ) break;
            }
         }
         else
         {
            for ( size_t l1 = 0; l1 < selectedL1Jets.size() ; ++l1 )
            {
               if ( jet -> deltaR(*selectedL1Jets[l1]) < matchonoffdrmax_ ) l1match = true;
               if ( l1match ) break;
            }
         }
         if ( !l1match && matchonoff_ ) continue;
         
         bool l2match = false;
         for ( size_t l2 = 0; l2 < selectedL2Jets.size() ; ++l2 )
         {
            if ( jet -> deltaR(*selectedL2Jets[l2]) < matchonoffdrmax_ ) l2match = true;
            if ( l2match ) break;
         }
         if ( !l2match && matchonoff_ ) continue;
         
         bool l3match = false;
         for ( size_t l3 = 0; l3 < selectedL3Jets.size() ; ++l3 )
         {
            if ( jet -> deltaR(*selectedL3Jets[l3]) < matchonoffdrmax_ ) l3match = true;
            if ( l3match ) break;
         }
         if ( !l3match && matchonoff_ ) continue;
      }
      
      // fill histograms
      h1["pt" ]    -> Fill(jet->pt() , psw);
      h1["eta"]    -> Fill(jet->eta(), psw);
      h1["phi"]    -> Fill(jet->phi(), psw);
      
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
} //end main

std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx)
{
   std::vector<TriggerObject *> seltos;
   
   // select trigger objects
   for ( int j = 0; j < tos->size() ; ++j )
   {
      TriggerObject * to = &tos->at(j);
      if ( ntomin_[idx] < 1 )
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
