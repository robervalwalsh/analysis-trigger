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

std::vector<Jet *> SelectedIdJets(std::shared_ptr< Collection<Jet> >);
void CreateHistograms(std::map<std::string, TH1F*> &);
void CreateHistograms(std::map<std::string, TH2F*> &);
template <typename T>
void WriteHistograms(std::map<std::string, T*> &);
bool JetJetCorrelations(const std::vector<Jet *> &);
bool JetOffOnlMatch(const Jet &, const std::string & exception = "");
float GetBTag(const Jet & jet, const std::string & algo);

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
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
      auto l = &obj - &triggerObjects_[0];
      if ( njetsmin_ == 1 )
      {
         triggerObjectsMatches_[l] = 1;
         triggerObjectsMatchesRank_[l] = 1;
      }
   }

   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   std::map<std::string, TH2F*> h2;
   CreateHistograms(h1);
   CreateHistograms(h2);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   
   // MAIN LOOP
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
// resets
      bool goodEvent = true;
      
      std::map<std::string,bool> isGood;
      for ( int j = 0; j < njetsmin_ ; ++j )
      {
         isGood[Form("pt%d",j)]    = true; 
         isGood[Form("eta%d",j)]   = true; 
         isGood[Form("btag%d",j)]  = true; 
         for ( int k = j+1; k < njetsmin_ && j < njetsmin_; ++k )
         {
            isGood[Form("dr%d%d",j,k)]   = true;
            isGood[Form("deta%d%d",j,k)] = true;
         }
      }
// -- resets
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
// read event
      analysis.event(i);
      
      // select run ranges
      if ( runmin_ > 0 && analysis.run() < runmin_ ) continue;
      if ( runmax_ > 0 && analysis.run() > runmax_ ) continue;
      
      // data-only
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
// Online
      // trigger fired
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
      // emulation of trigger should be called here
      
      // -- trigger fired
            
      // match offline to online - do not use in case of trigger emulation,
      // need to use jet->matchTo(trigobj, name, deltaR), if matched retrieve using jet->matched(name)
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects_,0.3);
// -- online
      
// Offline
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets = SelectedIdJets(slimmedJets);
      
      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
      // kinematic and offline btagging selection - n leading jets
      // *** DOING THIS WAY ONE HAS MORE CONTROL OF THE CONDITIONS OF THE SELECTION ***
      for ( int j = 0 ; j < njetsmin_ ; ++j )
      {
         if ( selectedJets[j] -> pt() < jetsptmin_[j]                       && !(jetsptmin_[j] < 0) )   isGood[Form("pt%d",j)]       = false;
         if ( fabs(selectedJets[j] -> eta()) > jetsetamax_[j]               && !(jetsetamax_[j] < 0) )  isGood[Form("eta%d",j)]      = false;
         if ( GetBTag(*selectedJets[j],btagalgo_) < jetsbtagmin_[j]         && !(jetsbtagmin_[j] < 0) ) isGood[Form("btag%d",j)]     = false;
         for ( int k = j+1; k < njetsmin_ && j < njetsmin_; ++k )
         {
            if ( selectedJets[j]->deltaR(*selectedJets[k]) < drmin_ )                                    isGood[Form("dr%d%d",j,k)]   = false;
         }
      }
      if ( njetsmin_ > 1 )
      {
         if ( fabs(selectedJets[0]->eta() - selectedJets[1]->eta()) > detamax_ && !(detamax_ < 0) )        isGood[Form("deta%d%d",0,1)] = false;
      }
      
      for ( auto & ok : isGood )
      {
         goodEvent = ( goodEvent && ok.second );
      }
      if ( ! goodEvent ) continue;
      
      // Trigger matching
      for ( size_t l = 0 ; l < triggerObjects_.size(); ++l )
      {
         if ( triggerObjectsMatches_[l] == triggerObjectsMatchesRank_[l] )
         {
            goodEvent = true;
            for ( int m = 0; m < triggerObjectsMatchesRank_[l]; ++m )
            {
               if ( ! selectedJets[m]->matched(triggerObjects_[l]) )
               {
                  goodEvent = false;
                  break;
               }
            }
            if ( ! goodEvent ) break;
         }
         else if ( triggerObjectsMatchesRank_[l] > 0 )
         {
            goodEvent = true;
            int m = triggerObjectsMatchesRank_[l]-1;
            if ( ! selectedJets[m]->matched(triggerObjects_[l]) )
            {
               goodEvent = false;
               break;
            }
         }
         else if ( triggerObjectsMatchesRank_[l] < 0 )
         {
            goodEvent = false;
            for ( int m = 0; m < abs(triggerObjectsMatchesRank_[l]) ; ++m )
            {
               goodEvent = (goodEvent || selectedJets[m]->matched(triggerObjects_[l]) );
            }
         }
      }
      
      if ( ! goodEvent ) continue;
      
      
      // Get prescales
      float psweight = 1.;
      if ( psweight_ && !isMC_ )
      {
         int hltps = analysis.triggerPrescale(hltPath_);
         int l1ps  = analysis.triggerPrescale(l1Seed_);
         psweight = (float)hltps*l1ps/psnorm_;
      }

      // Filling histograms...
      
      
      for ( int j = 0 ; j < njetsmin_ ; ++j )
      {
          h1[Form("pt_jet%d",j)]      -> Fill(selectedJets[j]->pt(),psweight);
          h1[Form("eta_jet%d",j)]     -> Fill(selectedJets[j]->eta(),psweight);
          h1[Form("btag_jet%d",j)]    -> Fill(GetBTag(*selectedJets[j],btagalgo_),psweight);
          for ( int k = j+1; k < njetsmin_ && j < njetsmin_; ++k )
          {
             h1[Form("dr_jet%d%d",j,k)] -> Fill(selectedJets[j]->deltaR(*selectedJets[k]),psweight);
             h1[Form("deta_jet%d%d",j,k)] -> Fill(fabs(selectedJets[j]->eta() - selectedJets[k]->eta()),psweight);
          }
      }
      
   }
   
   // Write histograms to output
   WriteHistograms(h1);
   WriteHistograms(h2);
   
   hout.Close();
   
// PRINT OUTS   
   std::cout << "AnalysisDoubleJet: program finished" << std::endl;
   
} //end main

// =====================================================================================
// =====================================================================================
// functions
// =====================================================================================
// =====================================================================================
std::vector<Jet *> SelectedIdJets(std::shared_ptr< Collection<Jet> > jets)
{
   std::vector<Jet *> seljets;
   for ( int j = 0 ; j < jets->size() ; ++j )
   {
      if ( jetsid_ == "TIGHT" ) // LOOSE is default if TIGHT not given
      {
         if ( jets->at(j).idTight() ) seljets.push_back(&jets->at(j));
      }
      else
      {
         if ( jets->at(j).idLoose() ) seljets.push_back(&jets->at(j));
      }
   }
   return seljets;
}

bool JetOffOnlMatch(const Jet & jet, const std::string & exception)
{
   bool  goodEvent = true;
   bool matched[10] = {true,true,true,true,true,true,true,true,true,true};
   for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
   {      
      if ( triggerObjects_[io] == exception ) continue;
      if ( ! jet.matched(triggerObjects_[io]) ) matched[io] = false;
      goodEvent = ( goodEvent && matched[io] );
   }
   return goodEvent;

}

void CreateHistograms(std::map<std::string, TH1F*> & h1)
{
   for ( int j = 0; j < njetsmin_; ++j )
   {
      h1[Form("pt_jet%d",j)]       = new TH1F(Form("pt_jet%d",j)     , "" ,100 , 0   , 1000  );
      h1[Form("eta_jet%d",j)]      = new TH1F(Form("eta_jet%d",j)    , "" , 60 , -3, 3 );
      h1[Form("btag_jet%d",j)]     = new TH1F(Form("btag_jet%d",j)   , "" , 100 , 0, 1 );
      for ( int k = j+1; k < njetsmin_ && j < njetsmin_; ++k )
      {
         h1[Form("dr_jet%d%d",j,k)]     = new TH1F(Form("dr_jet%d%d",j,k)     , "" , 50 , 0, 5 );
         h1[Form("deta_jet%d%d",j,k)]   = new TH1F(Form("deta_jet%d%d",j,k)   , "" ,100 , 0,10 );
      }
   }
   
}
void CreateHistograms(std::map<std::string, TH2F*> & h2)
{
   
}


template <typename T>
void WriteHistograms(std::map<std::string, T*> & h)
{
   // Write histograms to output
   for (auto & ih : h)
   {
      ih.second -> Write();
   }

}

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
   else if ( btagalgo_ == "deepbvsall" )
   {
      btag = jet.btag("btag_deepbvsall");
   }
   else
   {
      btag = -9999;
   }
   return btag;
}
