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

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Muons
   analysis.addTree<Muon> ("Muons","MssmHbb/Events/slimmedMuons");
   
   // Trigger path info
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   // Trigger objects
   for ( auto & obj : triggerObjects_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));
   }

   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   
   h1["n_muon"]   = new TH1F("n_muon"   , "" ,   5,  0   ,   5  );
   h1["pt_muon"]  = new TH1F("pt_muon"  , "" , 100,  0   ,  50  );
   h1["eta_muon"] = new TH1F("eta_muon" , "" ,  12, -2.4 ,   2.4);
   h1["phi_muon"] = new TH1F("phi_muon" , "" ,  16, -3.2 ,   3.2);
      
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
      
      // Muons
      auto slimmedMuons = analysis.collection<Muon>("Muons");
      // selection of muons
      std::vector<Muon *> selectedLooseMuons;
      for ( int m = 0 ; m < slimmedMuons->size() ; ++m )
      {
         Muon * muon = &(slimmedMuons->at(m));
         if ( muon -> isLooseMuon() )
            selectedLooseMuons.push_back(muon);
      }
      if ( (int)selectedLooseMuons.size() < nmuonsmin_ ) continue;
      std::vector<Muon *> selectedMuons;
      // selection only on the leading loose id'ed muons
      for ( int m = 0 ; m < (int)selectedLooseMuons.size() ; ++m )
      {
         Muon * muon = selectedLooseMuons.at(m);
         if ( muon->pt() > muonsptmin_[0] && fabs(muon->eta()) < muonsetamax_[0] )
            selectedMuons.push_back(muon);
      }
      if ( (int)selectedMuons.size() < nmuonsmin_ ) continue;
      
      Muon * muon = selectedMuons.at(0);

      
// Trigger analysis
// ==================================      
      bool hltAccept = analysis.triggerResult(hltPath_);
      bool l1Accept  = analysis.triggerResult(l1Seed_);
      if ( !hltAccept || !l1Accept ) continue; // Selecting one of the OR L1 seeds.
      
      std::vector<TriggerObject *> selectedL1Muons;
      std::vector<TriggerObject *> selectedL3Muons;
      // select L1 objects
      if ( triggerObjects_.size() > 1 )
      {
         auto l1muons = analysis.collection<TriggerObject>(triggerObjects_[0]);
         for ( int m = 0; m < l1muons->size() ; ++m )
         {
            TriggerObject * l1muon = &l1muons->at(m);
            if ( nl1muonsmin_ < 1 )
            {
               selectedL1Muons.push_back(l1muon);
            }
            else  // L1 emulation
            {
               if ( l1muon -> pt() >= l1muonsptmin_[0] && fabs( l1muon -> eta() ) <= l1muonsetamax_[0] ) 
                  selectedL1Muons.push_back(l1muon);
            }
         }
         if ( (int)selectedL1Muons.size() < nl1muonsmin_  ) continue;
      
         // select L3 objects
         auto l3muons = analysis.collection<TriggerObject>(triggerObjects_[1]);
         for ( int m = 0; m < l3muons->size() ; ++m )
         {
            TriggerObject * l3muon = &l3muons->at(m);
            if ( nl3muonsmin_ < 1 )
            {
               selectedL3Muons.push_back(l3muon);
            }
            else  // L3 emulation
            {
               if ( l3muon -> pt() >= l3muonsptmin_[0] && fabs( l3muon -> eta() ) <= l3muonsetamax_[0] ) 
                  selectedL3Muons.push_back(l3muon);
            }
         }
         if ( (int)selectedL3Muons.size() < nl3muonsmin_  ) continue;
      }
      
      // match leading muon to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
      if ( triggerObjects_.size() > 1 )
      {
         bool l1match = false;
         for ( size_t l1 = 0; l1 < selectedL1Muons.size() ; ++l1 )
         {
            if ( muon -> deltaR(*selectedL1Muons[l1]) < matchonoffdrmax_ ) l1match = true;
            if ( l1match ) break;
         }
         if ( !l1match && matchonoff_ ) continue;
         
         bool l3match = false;
         for ( size_t l3 = 0; l3 < selectedL3Muons.size() ; ++l3 )
         {
            if ( muon -> deltaR(*selectedL3Muons[l3]) < matchonoffdrmax_ ) l3match = true;
            if ( l3match ) break;
         }
         if ( !l3match && matchonoff_ ) continue;
      }
      
      // fill histograms
      h1["pt_muon" ]    -> Fill(muon->pt() );
      h1["eta_muon"]    -> Fill(muon->eta());
      h1["phi_muon"]    -> Fill(muon->phi());
      
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
} //end main

