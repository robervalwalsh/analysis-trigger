/*

Thhis macro plots offline leading and second leading jet information from a double jet trigger

*/

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

float GetBTag(const Jet & jet, const std::string & algo);

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // some specific requirements
   if ( njetsmin_ != 2 )
   {
      std::cout << "*** This is a double jet analysis, nJetsMin must be equal to two. ***" << std::endl;
      return -1;
   }
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Jets
   analysis.addTree<Jet> ("Jets",jetsCol_);
   // Trigger path info
   analysis.triggerResults(triggerCol_);
   // Trigger objects
   for ( auto & obj : triggerObjects_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("%s/%s",triggerObjDir_.c_str(),obj.c_str()));
   }

   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   std::map<std::string, TH2F*> h2;
   int nbins = 26;
   double ptbins[27] = {20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,350,400,500,700,1000};
   
   h1["pt_jet0"  ]      = new TH1F("pt_jet0"      , "" , 100, 0   , 1000);
   h1["pt_jet0_var"]    = new TH1F("pt_jet0_var"  , "" , nbins, ptbins);
   h1["eta_jet0" ]      = new TH1F("eta_jet0"     , "" , 11 , -2.2, 2.2 );
   h1["phi_jet0" ]      = new TH1F("phi_jet0"     , "" , 16,  -3.2, 3.2 );
   h1["btag_jet0"]      = new TH1F("btag_jet0"    , "" , 200, 0   , 1   );
   h1["btaglog_jet0"]   = new TH1F("btaglog_jet0" , "" , 20 , 0   , 10  );
   
   h1["pt_jet1"  ]    = new TH1F("pt_jet1"      , "" , 100, 0   , 1000);
   h1["pt_jet1_var"]  = new TH1F("pt_jet1_var"  , "" , nbins, ptbins);
   h1["eta_jet1" ]    = new TH1F("eta_jet1"     , "" , 11 , -2.2, 2.2 );
   h1["phi_jet1" ]    = new TH1F("phi_jet1"     , "" , 16,  -3.2, 3.2 );
   h1["btag_jet1"]    = new TH1F("btag_jet1"    , "" , 200, 0   , 1   );
   h1["btaglog_jet1"] = new TH1F("btaglog_jet1" , "" , 20 , 0   , 10  );
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   int nsel[10] = { };

   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      bool goodEvent = true;
      
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
      
      ++nsel[0];
      
      // match offline to online
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects_,0.3);
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         if ( jetsid_ == "TIGHT" ) // LOOSE is default if TIGHT not given
         {
            if ( slimmedJets->at(j).idTight() ) selectedJets.push_back(&slimmedJets->at(j));
         }
         else
         {
            if ( slimmedJets->at(j).idLoose() ) selectedJets.push_back(&slimmedJets->at(j));
         }
      }
      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
      ++nsel[1];
      
      // kinematic and offline btagging selection - n leading jets
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         if ( jet->pt() < jetsptmin_[j] || fabs(jet->eta()) > jetsetamax_[j] )
         {
            goodEvent = false;
            break;
         }
      }
      if ( ! goodEvent ) continue;
      
      // deltaR
      for ( int j1 = 0; j1 < njetsmin_-1; ++j1 )
      {
         const Jet & jet1 = *selectedJets[j1];
         for ( int j2 = j1+1; j2 < njetsmin_; ++j2 )
         {
            const Jet & jet2 = *selectedJets[j2];
            if ( fabs(jet1.deltaR(jet2)) < drmin_ ) goodEvent = false;
         }
      }
      
      if ( ! goodEvent ) continue;
      
      // jet0 must be matched
      bool jet0matched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
      {       
         if ( ! selectedJets[1]->matched(triggerObjects_[io]) ) jet0matched[io] = false;
         goodEvent = ( goodEvent && jet0matched[io] );
      }
      if ( ! goodEvent ) continue;
      
      // is jet1 jet matched?
      bool jet1matched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
      {       
         if ( ! selectedJets[0]->matched(triggerObjects_[io]) ) jet1matched[io] = false;
         goodEvent = ( goodEvent && jet1matched[io] );
      }
      
      if ( ! goodEvent ) continue;
      
      ++nsel[2];
     
      // fill histograms for tag jet 
      
      // fill histograms for probe jet not matched yet 
      h1["pt_jet0"  ]    -> Fill(selectedJets[0]->pt()  );
      h1["pt_jet0_var"]  -> Fill(selectedJets[0]->pt()  );
      h1["eta_jet0" ]    -> Fill(selectedJets[0]->eta() );
      h1["phi_jet0" ]    -> Fill(selectedJets[0]->phi() );
      h1["btag_jet0"]    -> Fill(GetBTag(*selectedJets[0],btagalgo_));
      h1["btaglog_jet0"] -> Fill(-log(1-GetBTag(*selectedJets[0],btagalgo_)));
      
      h1["pt_jet1"  ]    -> Fill(selectedJets[1]->pt()  );
      h1["pt_jet1_var"]  -> Fill(selectedJets[1]->pt()  );
      h1["eta_jet1" ]    -> Fill(selectedJets[1]->eta() );
      h1["phi_jet1" ]    -> Fill(selectedJets[1]->phi() );
      h1["btag_jet1"]    -> Fill(GetBTag(*selectedJets[1],btagalgo_));
      h1["btaglog_jet1"] -> Fill(-log(1-GetBTag(*selectedJets[1],btagalgo_)));
      
      
   }
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   for (auto & ih2 : h2)
   {
      ih2.second -> Write();
   }
   
   
   hout.Close();
   
// PRINT OUTS   
   
   std::cout << "SingleJetTriggerAnalysis: program finished" << std::endl;
   
   
      
   
} //end main

float GetBTag(const Jet & jet, const std::string & algo)
{
   float btag;
   if ( btagalgo_ == "csvivf" || btagalgo_ == "csv" )
   {
      btag = jet.btag("btag_csvivf");
   }
   else if ( btagalgo_ == "deepcsv" )
   {
      btag = jet.btag("btag_deepb") + jet.btag("btag_deepbb");
   }
   else
   {
      btag = -9999;
   }
   return btag;
}
