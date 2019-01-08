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
bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos);

bool IsGoodPrimaryVertex(const Vertex & pv);

bool TriggerAccept(Analysis & analysis, const Muon & muon, const std::string & type = "nom");
bool TriggerAccept(Analysis & analysis, const Jet & jet, const std::string & type = "nom");


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
   
   // HLT_Mu8
//   int nbins = 18;
//   double ptbins[19] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,30,50,100};

   // HLT_Mu12
   float nbins = 20;
   float ptbins[21] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,35,60,120};
      
   // numerator histograms
   h1["n_num"]       = new TH1F("n_num"   , "" ,   5,  0   ,   5  );
   h1["pt_num"]      = new TH1F("pt_num"  , "" , 150,  0   ,  150 );
   h1["pt_var_num"]  = new TH1F("pt_var_num"  , "" , nbins, ptbins );
   h1["eta_num"] = new TH1F("eta_num" , "" , 100, -2.5 ,   2.5);
   h1["phi_num"] = new TH1F("phi_num" , "" , 128, -3.2 ,   3.2);
   
   h1["n_den"]       = new TH1F("n_den"   , "" ,   5,  0   ,   5  );
   h1["pt_den"]      = new TH1F("pt_den"  , "" , 150,  0   ,  150 );
   h1["pt_var_den"]  = new TH1F("pt_var_den"  , "" , nbins, ptbins );
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
      
      Muon * muon = selectedMuons.at(0);
//      std::cout << i<< std::endl;
//      for ( auto & mu: selectedMuons )
//         std::cout << " offline: " << mu->pt() << "  " << mu-> eta() << "   " << mu->phi() << "  " << mu->isMediumMuon() << "  "  << std::endl;
//      std::cout << "----------" << std::endl;
      
      
// Trigger analysis, reference trigger (denominator)
// =================================================  
      bool refAccept = TriggerAccept(analysis,*muon,"ref");
      if ( ! refAccept ) continue;
      
      float pswref = 1;
//      if ( psweight_ ) pswref = analysis.triggerPrescale(hltPathRef_)*analysis.triggerPrescale(l1SeedRef_);
      
      // fill histograms
      h1["pt_den"]      -> Fill(muon->pt() ,pswref);
      h1["pt_var_den"]  -> Fill(muon->pt() ,pswref);
      h1["eta_den"]     -> Fill(muon->eta(),pswref);
      h1["phi_den"]     -> Fill(muon->phi(),pswref);
      
// // Trigger analysis, nominal trigger (numerator)
// // =================================================    
      bool nomAccept = TriggerAccept(analysis,*muon,"nom");
      
      if ( ! nomAccept ) continue;
      
      float psw = 1;
//      if ( psweight_ ) psw = analysis.triggerPrescale(hltPath_)*analysis.triggerPrescale(l1Seed_);
      if ( psweight_ ) psw = analysis.triggerPrescale(hltPath_);
      
      // fill histograms
      h1["pt_num"]      -> Fill(muon->pt() , psw);
      h1["pt_var_num"]  -> Fill(muon->pt() , psw);
      h1["eta_num"]     -> Fill(muon->eta(), psw);
      h1["phi_num"]     -> Fill(muon->phi(), psw);
      
//      std::cout << "========" << std::endl;
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
   std::cout << "SingleMuonTriggerEfficiencies: program finished" << std::endl;
   
} //end main


// ================================

bool TriggerAccept(Analysis & analysis, const Muon & muon, const std::string & type)
{
   
   std::string hlt = hltPath_;
   std::string l1  = l1Seed_;
   int l1tnmin = l1tmuonsnmin_;
   std::vector<float> l1tptmin = l1tmuonsptmin_;
   std::vector<float> l1tetamax = l1tmuonsetamax_;
   
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
      l1tnmin = l1tmuonsrefnmin_;
      l1tptmin = l1tmuonsrefptmin_;
      l1tetamax = l1tmuonsrefetamax_;
      tonmin = torefnmin_;
      toptmin = torefptmin_;
      toetamax = torefetamax_;
      tos = triggerObjectsRef_;
   }
   
// Trigger analysis, reference trigger (denominator)
// =================================================      
   bool hltAccept = analysis.triggerResult(hlt);
   // see SingleMuonTriggerAnalysis.cc for the usage of prescales if that is needed
   
   bool l1Accept  = true;
   if ( l1 != "" )  l1Accept = analysis.triggerResult(l1);
   
//   std::cout << type << "   "  << hltAccept << "  " << l1Accept << std::endl;
   
   if ( !hltAccept || !l1Accept ) return false; // Selecting one of the OR L1 seeds.
   
   // Emulation with L1T candidate
   // L1TMuons
   // selection of muons
   std::vector<L1TMuon *> selectedL1TMuons;
   if ( l1tnmin > 0 )
   {
      auto l1tmuons = analysis.collection<L1TMuon>("L1TMuons");
      for ( int j = 0 ; j < l1tmuons->size() ; ++j )
      {
         L1TMuon * l1tmuon = &(l1tmuons->at(j));
         if ( l1tmuon->pt() > l1tptmin[0] && fabs(l1tmuon->eta()) < l1tetamax[0] )
            selectedL1TMuons.push_back(l1tmuon);
      }
      if ( (int)selectedL1TMuons.size() < l1tnmin ) return false;
   }
   
   // Emulation with trigger objects
   std::vector<TriggerObject *> selectedL1Muons;
   std::vector<TriggerObject *> selectedL2Muons;
   std::vector<TriggerObject *> selectedL3Muons;
   
   if ( tos.size() > 2 )
   {
      // select L1 objects
      auto l1muons = analysis.collection<TriggerObject>(tos[0]);
      selectedL1Muons = SelectTriggerObjects(l1muons,0,tonmin,toptmin,toetamax);
//      if ( (int)selectedL1Muons.size() < tonmin[0]  ) return false;
      
      // select L2 objects
      auto l2muons = analysis.collection<TriggerObject>(tos[1]);
      selectedL2Muons = SelectTriggerObjects(l2muons,1,tonmin,toptmin,toetamax);
//      if ( (int)selectedL2Muons.size() < tonmin[1]  ) return false;
   
      // select L3 objects
      auto l3muons = analysis.collection<TriggerObject>(tos[2]);
      selectedL3Muons = SelectTriggerObjects(l3muons,2,tonmin,toptmin,toetamax);
//      if ( (int)selectedL3Muons.size() < tonmin[2]  ) return false;
   }
   
//    if ( type != "ref" )
//    {
//       for ( auto & ol1 : selectedL1Muons )
//          std::cout << " L1 "  << ol1->pt() << "  " << ol1->eta() << "  " << ol1->phi() << std::endl;
//       for ( auto & ol2 : selectedL2Muons )
//          std::cout << " L2 "  << ol2->pt() << "  " << ol2->eta() << "  " << ol2->phi() << std::endl;
//       for ( auto & ol3 : selectedL3Muons )
//          std::cout << " L3 "  << ol3->pt() << "  " << ol3->eta() << "  " << ol3->phi() << std::endl;
//       
//    }
      if ( (int)selectedL1Muons.size() < tonmin[0]  ) return false;
      if ( (int)selectedL2Muons.size() < tonmin[1]  ) return false;
      if ( (int)selectedL3Muons.size() < tonmin[2]  ) return false;
      
//       if ( type == "nom" ) std::cout << "trigger fired" << std::endl;

   


   // match leading muon to trigger objects - will be done by hand to be sure to use the correct emulated trigger objects
   if ( matchonoff_ )
   {
      bool l1match = false;
      if ( tonmin[0] > 0 )       l1match = MatchOnlineOffline(muon,selectedL1Muons);
      else if (l1tnmin > 0)      l1match = MatchOnlineOffline(muon,selectedL1TMuons);
      else l1match = true;
      if ( !l1match ) return false;
      
      bool l2match = MatchOnlineOffline(muon,selectedL2Muons);
      if ( !l2match ) return false;
      
      bool l3match = MatchOnlineOffline(muon,selectedL3Muons);
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


bool IsGoodPrimaryVertex(const Vertex & pv)
{
   bool good = false;
   
   if ( ! pv.fake() && pv.ndof() > 4 && fabs(pv.z()) <= 24 && pv.rho() <=2 ) good = true;
   else good = false;
   
   return good;
}
