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
void CreateHistograms(std::map<std::string, TH1F*> &, std::map<std::string, TH2F*> &);
void WriteHistograms(std::map<std::string, TH1F*> &, std::map<std::string, TH2F*> &);
bool JetJetCorrelations(const std::vector<Jet *> &);
bool JetOffOnlMatch(const Jet &, const std::string & exception = "");
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
      std::cout << "*** This is a tag and probe between two jets, nJetsMin must be equal to two. ***" << std::endl;
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
   CreateHistograms(h1,h2);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
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
      
      // match offline to online
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects_,0.5);
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets = SelectedIdJets(slimmedJets);
      
      
      if ( (int)selectedJets.size() < njetsmin_ || (int)selectedJets.size() > njetsmax_ ) continue;
//      if ( (int)selectedJets.size() < njetsmin_  ) continue;
      
      // additional jet veto
      if ( (int)selectedJets.size() > njetsmin_ && ptmax_.size() > 0 )
      {
         if ( selectedJets[njetsmin_]->pt() > ptmax_[0] ) continue;
      }
      
      // kinematic and offline btagging selection - n leading jets
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         if ( jet->pt() < jetsptmin_[j] || jet->pt() > jetsptmax_[j] || fabs(jet->eta()) > jetsetamax_[j] || GetBTag(*jet,btagalgo_) < jetsbtagmin_[j] )
         {
            goodEvent = false;
            break;
         }
      }
      if ( ! goodEvent ) continue;
      
      // deltaR and deltaPhi and imbalance
      if ( ! JetJetCorrelations(selectedJets) ) continue;
      // tag jets must be matched
      if ( ! JetOffOnlMatch(*(selectedJets[1])) ) continue;
      // probe matched except for the btag trigger
      if ( ! JetOffOnlMatch(*(selectedJets[0]),triggerObjects_.back()) ) continue;
      
      // Filling histograms...
      h1["n_jets"] -> Fill(selectedJets.size());
      if ( (int)selectedJets.size() > njetsmin_ )
      {
         h1["pt_jet3"] -> Fill( selectedJets[njetsmin_]->pt() );
      }
      h1["delta_pt"] -> Fill((selectedJets[0]->pt()-selectedJets[1]->pt())/selectedJets[0]->pt());
      h1["delta_phi"] -> Fill(fabs(selectedJets[0]->deltaPhi(*selectedJets[1])));
      
      // fill histograms for tag jet 
      h1["pt_tag"  ]    -> Fill(selectedJets[1]->pt()  );
      h1["pt_tag_var"]  -> Fill(selectedJets[1]->pt()  );
      h1["eta_tag" ]    -> Fill(selectedJets[1]->eta() );
      h1["phi_tag" ]    -> Fill(selectedJets[1]->phi() );
      h1["btag_tag"]    -> Fill(GetBTag(*selectedJets[1],btagalgo_));
      h1["btaglog_tag"] -> Fill(-log(1-GetBTag(*selectedJets[1],btagalgo_)));
      
      // fill histograms for probe jet not matched yet 
      h1["pt_probe_den"  ]    -> Fill(selectedJets[0]->pt()  );
      h1["pt_probe_den_var"]  -> Fill(selectedJets[0]->pt()  );
      h1["eta_probe_den" ]    -> Fill(selectedJets[0]->eta() );
      h1["phi_probe_den" ]    -> Fill(selectedJets[0]->phi() );
      h1["btag_probe_den"]    -> Fill(GetBTag(*selectedJets[0],btagalgo_));
      h1["btaglog_probe_den"] -> Fill(-log(1-GetBTag(*selectedJets[0],btagalgo_)));
      
      h1["qglikelihood_tag"] -> Fill(selectedJets[1]->qgLikelihood());
      h1["qglikelihood_probe_den"] -> Fill(selectedJets[0]->qgLikelihood());
      h2["qglikelihood_tag_x_probe_den"] -> Fill(selectedJets[1]->qgLikelihood(),selectedJets[0]->qgLikelihood()); 
      
      if ( selectedJets[1]->qgLikelihood() > 0.8 )
      {
         h1["qglikelihood_tag0p8"] -> Fill(selectedJets[1]->qgLikelihood());
         h1["qglikelihood_probe_den_tag0p8"] -> Fill(selectedJets[0]->qgLikelihood());
      }
       
      if ( selectedJets[0]->qgLikelihood() > 0.8 )
      {
         h1["qglikelihood_tag_probe0p8"] -> Fill(selectedJets[1]->qgLikelihood());
         h1["qglikelihood_probe_den_0p8"] -> Fill(selectedJets[0]->qgLikelihood());
      }
      
      // probe jets must be matched now
      if ( ! JetOffOnlMatch(*(selectedJets[0])) ) continue;
      
      // fill histograms for probe jet 
      h1["pt_probe_num"  ]    -> Fill(selectedJets[0]->pt()  );
      h1["pt_probe_num_var"]  -> Fill(selectedJets[0]->pt()  );
      h1["eta_probe_num" ]    -> Fill(selectedJets[0]->eta() );
      h1["phi_probe_num" ]    -> Fill(selectedJets[0]->phi() );
      h1["btag_probe_num"]    -> Fill(GetBTag(*selectedJets[0],btagalgo_));
      h1["btaglog_probe_num"] -> Fill(-log(1-GetBTag(*selectedJets[0],btagalgo_)));
      h1["qglikelihood_probe_num"] -> Fill(selectedJets[0]->qgLikelihood());
      
   }
   
   // Write histograms to output
   WriteHistograms(h1,h2);
   hout.Close();
   
// PRINT OUTS   
   std::cout << "SingleJetTriggerAnalysis: program finished" << std::endl;
   
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

bool JetJetCorrelations(const std::vector<Jet *> & selectedJets)
{
   bool goodEvent = true;
   for ( int j1 = 0; j1 < njetsmin_-1; ++j1 )
   {
      const Jet & jet1 = *selectedJets[j1];
      for ( int j2 = j1+1; j2 < njetsmin_; ++j2 )
      {
         const Jet & jet2 = *selectedJets[j2];
         if ( fabs(jet1.deltaR(jet2)) < drmin_ ) goodEvent = false;
         if ( fabs(jet1.deltaPhi(jet2)) < dphimin_ ) goodEvent = false;
         if ( fabs(jet1.pt()-jet2.pt())/jet1.pt() > ptimbalmax_ ) goodEvent = false;
      }
   }
   return goodEvent;

}

void CreateHistograms(std::map<std::string, TH1F*> & h1, std::map<std::string, TH2F*> & h2)
{
//   int nbins = 12+7+6+2;
//   double ptbins[28] = {40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,350,400,450,500,550,600,800,1000};
   int nbins = 12+7+5+6;
//   double ptbins[31] = {10,15,20,25,30,35,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,350,400,500,700,1000};
   double ptbins[31] = {10,15,20,25,30,35,40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,320,380,450,550,700,1000};
   
   h1["n_jets"]        = new TH1F("n_jets"      , "" , 20 , 0   , 20  );
   h1["pt_jet3"]       = new TH1F("pt_jet3"     , "" ,100 , 0   , 1000  );
   
   h1["pt_tag"  ]      = new TH1F("pt_tag"      , "" , 100, 0   , 1000);
   h1["pt_tag_var"]    = new TH1F("pt_tag_var"  , "" , nbins, ptbins);
   h1["eta_tag" ]      = new TH1F("eta_tag"     , "" , 11 , -2.2, 2.2 );
   h1["phi_tag" ]      = new TH1F("phi_tag"     , "" , 16,  -3.2, 3.2 );
   h1["btag_tag"]      = new TH1F("btag_tag"    , "" , 200, 0   , 1   );
   h1["btaglog_tag"]   = new TH1F("btaglog_tag" , "" , 20 , 0   , 10  );
   
   h1["qglikelihood_tag"] = new TH1F("qglikelihood_tag","",100, 0,1);
   h1["qglikelihood_tag0p8"] = new TH1F("qglikelihood_tag0p8","",100, 0,1);
   h1["qglikelihood_tag_probe0p8"] = new TH1F("qglikelihood_tag_probe0p8","",100, 0,1);
   
   h1["pt_probe_den"  ]    = new TH1F("pt_probe_den"      , "" , 100, 0   , 1000);
   h1["pt_probe_den_var"]  = new TH1F("pt_probe_den_var"  , "" , nbins, ptbins);
   h1["eta_probe_den" ]    = new TH1F("eta_probe_den"     , "" , 11 , -2.2, 2.2 );
   h1["phi_probe_den" ]    = new TH1F("phi_probe_den"     , "" , 16,  -3.2, 3.2 );
   h1["btag_probe_den"]    = new TH1F("btag_probe_den"    , "" , 200, 0   , 1   );
   h1["btaglog_probe_den"] = new TH1F("btaglog_probe_den" , "" , 20 , 0   , 10  );
   
   
   h1["delta_pt"] = new TH1F("delta_pt","", 100, 0, 1);
   h1["delta_phi"] = new TH1F("delta_phi","", 320, 0, 3.2);

   h1["qglikelihood_probe_den"] = new TH1F("qglikelihood_probe_den","",100, 0,1);
   h1["qglikelihood_probe_den_tag0p8"] = new TH1F("qglikelihood_probe_den_tag0p8","",100, 0,1);
   h1["qglikelihood_probe_den_0p8"] = new TH1F("qglikelihood_probe_den_0p8","",100, 0,1);
   
   h1["pt_probe_num"  ]    = new TH1F("pt_probe_num"      , "" , 100, 0   , 1000);
   h1["pt_probe_num_var"]  = new TH1F("pt_probe_num_var"  , "" , nbins, ptbins);
   h1["eta_probe_num" ]    = new TH1F("eta_probe_num"     , "" , 11 , -2.2, 2.2 );
   h1["phi_probe_num" ]    = new TH1F("phi_probe_num"     , "" , 16,  -3.2, 3.2 );
   h1["btag_probe_num"]    = new TH1F("btag_probe_num"    , "" , 200, 0   , 1   );
   h1["btaglog_probe_num"] = new TH1F("btaglog_probe_num" , "" , 20 , 0   , 10  );
      
   h1["qglikelihood_probe_num"] = new TH1F("qglikelihood_probe_num","",100, 0,1);
   
   
   h2["qglikelihood_tag_x_probe_den"] = new TH2F("qglikelihood_tag_x_probe_den","",200, 0,1,200, 0,1);

}

void WriteHistograms(std::map<std::string, TH1F*> & h1, std::map<std::string, TH2F*> & h2)
{
   // Write histograms to output
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   for (auto & ih2 : h2)
   {
      ih2.second -> Write();
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
   else
   {
      btag = -9999;
   }
   return btag;
}
