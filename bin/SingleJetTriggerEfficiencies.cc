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
bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos);

bool TriggerAccept(Analysis & analysis, const Jet & jet, const std::string & type = "nom");

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
   // merge nominal and reference
   std::vector <std::string> toall;
   toall.reserve( triggerObjects_.size() + triggerObjectsRef_.size() ); // preallocate memory
   toall.insert( toall.end(), triggerObjects_.begin(), triggerObjects_.end() );
   toall.insert( toall.end(), triggerObjectsRef_.begin(), triggerObjectsRef_.end() );
   // remove duplicates
   sort( toall.begin(), toall.end() );
   toall.erase( unique( toall.begin(), toall.end() ), toall.end() );
   
   for ( auto & obj : toall )
   {
      analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   }
   
   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   
   // numerator histograms
   h1["n_num"]   = new TH1F("n_num"   , "" ,   5,  0   ,   5  );
   h1["pt_num"]  = new TH1F("pt_num"  , "" , 500,  0   ,  1000  );
   h1["eta_num"] = new TH1F("eta_num" , "" , 100, -2.5 ,   2.5);
   h1["phi_num"] = new TH1F("phi_num" , "" , 128, -3.2 ,   3.2);
   
   h1["n_den"]   = new TH1F("n_den"   , "" ,   5,  0   ,   5  );
   h1["pt_den"]  = new TH1F("pt_den"  , "" , 500,  0   ,  1000  );
   h1["eta_den"] = new TH1F("eta_den" , "" , 100, -2.5 ,   2.5);
   h1["phi_den"] = new TH1F("phi_den" , "" , 128, -3.2 ,   3.2);
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%10000==0 ) std::cout << i << "  events processed! " << std::endl;
      
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
      
      
// Trigger analysis, reference trigger (denominator)
// =================================================  
      bool refAccept = TriggerAccept(analysis,*jet,"ref");
      if ( ! refAccept ) continue;
      
      float pswref = 1;
//      if ( psweight_ ) pswref = analysis.triggerPrescale(hltPathRef_)*analysis.triggerPrescale(l1SeedRef_);
      
      // fill histograms
      h1["pt_den" ]     -> Fill(jet->pt() ,pswref);
      h1["eta_den"]     -> Fill(jet->eta(),pswref);
      h1["phi_den"]     -> Fill(jet->phi(),pswref);
      
// // Trigger analysis, nominal trigger (numerator)
// // =================================================    
      bool nomAccept = TriggerAccept(analysis,*jet,"nom");
      
      if ( ! nomAccept ) continue;
      
      float psw = 1;
//      if ( psweight_ ) psw = analysis.triggerPrescale(hltPath_)*analysis.triggerPrescale(l1Seed_);
      if ( psweight_ ) psw = analysis.triggerPrescale(hltPath_);
      
      // fill histograms
      h1["pt_num" ]     -> Fill(jet->pt() , psw);
      h1["eta_num"]     -> Fill(jet->eta(), psw);
      h1["phi_num"]     -> Fill(jet->phi(), psw);
      
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
   std::cout << "SingleJetTriggerEfficiencies: program finished" << std::endl;
   
} //end main


// ================================

bool TriggerAccept(Analysis & analysis, const Jet & jet, const std::string & type)
{
   
   std::string hlt = hltPath_;
   std::string l1  = l1Seed_;
   int l1tnmin = l1tjetsnmin_;
   std::vector<float> l1tptmin = l1tjetsptmin_;
   std::vector<float> l1tetamax = l1tjetsetamax_;
   
   std::vector<std::string> tos = triggerObjects_;
   
   int  * tonmin = tonmin_;
   std::vector<float> * toptmin = toptmin_;
   std::vector<float> * toetamax = toetamax_;
   
   
   if ( trigemul_ && type == "nom" )
   {
      // emulating on top of the reference trigger, kind o f redundant but will keep it
      hlt = hltPathRef_;
      l1 = l1SeedRef_;
      tos = triggerObjectsRef_;
   }
   if ( type == "ref" )
   {
      hlt = hltPathRef_;
      l1 = l1SeedRef_;
      l1tnmin = l1tjetsrefnmin_;
      l1tptmin = l1tjetsrefptmin_;
      l1tetamax = l1tjetsrefetamax_;
      tonmin = torefnmin_;
      toptmin = torefptmin_;
      toetamax = torefetamax_;
      tos = triggerObjectsRef_;
   }
   
// Trigger analysis, reference trigger (denominator)
// =================================================      
   bool hltAccept = analysis.triggerResult(hlt);
   // see SingleJetTriggerAnalysis.cc for the usage of prescales if that is needed
   
   bool l1Accept  = true;
   if ( l1 != "" )  analysis.triggerResult(l1);
   
   if ( !hltAccept || !l1Accept ) return false; // Selecting one of the OR L1 seeds.
   
   // Emulation with L1T candidate
   // L1TJets
   // selection of jets
   std::vector<L1TJet *> selectedL1TJets;
   if ( l1tnmin > 0 )
   {
      auto l1tjets = analysis.collection<L1TJet>("L1TJets");
      for ( int j = 0 ; j < l1tjets->size() ; ++j )
      {
         L1TJet * l1tjet = &(l1tjets->at(j));
         if ( l1tjet->pt() > l1tptmin[0] && fabs(l1tjet->eta()) < l1tetamax[0] )
            selectedL1TJets.push_back(l1tjet);
      }
      if ( (int)selectedL1TJets.size() < l1tnmin ) return false;
   }
   
   // Emulation with trigger objects
   std::vector<TriggerObject *> selectedL1Jets;
   std::vector<TriggerObject *> selectedL2Jets;
   std::vector<TriggerObject *> selectedL3Jets;
   
   if ( tos.size() > 2 )
   {
      // select L1 objects
      auto l1jets = analysis.collection<TriggerObject>(tos[0]);
      selectedL1Jets = SelectTriggerObjects(l1jets,0,tonmin,toptmin,toetamax);
      if ( (int)selectedL1Jets.size() < tonmin[0]  ) return false;
      
      // select L2 objects
      auto l2jets = analysis.collection<TriggerObject>(tos[1]);
      selectedL2Jets = SelectTriggerObjects(l2jets,1,tonmin,toptmin,toetamax);
      if ( (int)selectedL2Jets.size() < tonmin[1]  ) return false;
   
      // select L3 objects
      auto l3jets = analysis.collection<TriggerObject>(tos[2]);
      selectedL3Jets = SelectTriggerObjects(l3jets,2,tonmin,toptmin,toetamax);
      if ( (int)selectedL3Jets.size() < tonmin[2]  ) return false;
   }

   // match leading jet to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
   if ( matchonoff_ )
   {
      bool l1match = false;
      if ( tonmin[0] > 0 )       l1match = MatchOnlineOffline(jet,selectedL1Jets);
      else if (l1tnmin > 0)      l1match = MatchOnlineOffline(jet,selectedL1TJets);
      else l1match = true;
      if ( !l1match ) return false;
      
      bool l2match = MatchOnlineOffline(jet,selectedL2Jets);
      if ( !l2match ) return false;
      
      bool l3match = MatchOnlineOffline(jet,selectedL3Jets);
      if ( !l3match ) return false;
   }
   
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


