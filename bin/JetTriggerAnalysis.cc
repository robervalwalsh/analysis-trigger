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


std::vector<L1TJet *> selectedL1TObjects_;
std::vector< std::vector<TriggerObject *> > selectedTriggerObjects_;


void SelectTriggerObjects(Analysis & analysis);
std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx);
std::vector<L1TJet *>        SelectTriggerObjects(const std::shared_ptr< Collection<L1TJet> > & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<TriggerObject *> & tos);
bool MatchOnlineOffline(const Jet & jet, const std::vector<L1TJet *> & tos);

bool IsGoodPrimaryVertex(const Vertex & pv);
bool TriggerAccept(Analysis & analysis);

void CreateHistograms(std::map<std::string, TH1F*> & h, const std::string & label, const int & nmin = 1);
void FillHistograms(std::map<std::string, TH1F*> & h, const std::string & label, const std::vector<TriggerObject *> & tos, const int & nmin = 1);
void FillHistograms(std::map<std::string, TH1F*> & h, const std::string & label, const std::vector<L1TJet *> & tos       , const int & nmin = 1);

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
   // Trigger path info
   analysis.triggerResults(Form("%s/%s",treePath_.c_str(),triggerCol_.c_str()));
   // Trigger objects
   for ( auto & obj : triggerObjects_ )  analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   
   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   CreateHistograms(h1,"l1t",l1tjetsnmin_);
   CreateHistograms(h1,"l1" ,tonmin_[0]);
   CreateHistograms(h1,"l2" ,tonmin_[1]);
   CreateHistograms(h1,"l3" ,tonmin_[2]);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%10000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i); // read event
      
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      
// Trigger selection 
// =================
      if ( ! TriggerAccept(analysis) ) continue;
      
      FillHistograms(h1, "l1t", selectedL1TObjects_       , l1tjetsnmin_);
      FillHistograms(h1, "l1" , selectedTriggerObjects_[0], tonmin_[0]);
      FillHistograms(h1, "l2" , selectedTriggerObjects_[1], tonmin_[1]);
      FillHistograms(h1, "l3" , selectedTriggerObjects_[2], tonmin_[2]);
      
   } // end of event loop
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
   std::cout << "program finished" << std::endl;
   
} //end main


// ================================

bool TriggerAccept(Analysis & analysis)
{
   // trigger bits
   bool l1Accept   = true;
   bool hltAccept  = true;
   bool emulAccept = true;
   if ( l1Seed_  != "" ) l1Accept  = analysis.triggerResult(l1Seed_);
   if ( hltPath_ != "" ) hltAccept = analysis.triggerResult(hltPath_);
   
   SelectTriggerObjects(analysis);
   // trigger emulator
   if ( ! trigemul_ ) return (l1Accept && hltAccept);
   
   return (l1Accept && hltAccept && emulAccept);

}


void SelectTriggerObjects(Analysis & analysis)
{
   selectedL1TObjects_.clear();
   selectedTriggerObjects_.clear();
   
   auto l1t = analysis.collection<L1TJet>("L1TJets");
   selectedL1TObjects_ = SelectTriggerObjects(l1t);
   
   std::vector< std::shared_ptr< Collection<TriggerObject> > > triggerObjects;
   for ( size_t i = 0 ; i < triggerObjects_.size() ; ++i )
   {
      triggerObjects.push_back(analysis.collection<TriggerObject>(triggerObjects_.at(i)));
      selectedTriggerObjects_.push_back(SelectTriggerObjects(triggerObjects.back(),i));
   }
}

std::vector<TriggerObject *> SelectTriggerObjects(const std::shared_ptr< Collection<TriggerObject> > & tos, const int & idx)
{
   std::vector<TriggerObject *> seltos;
   
   // select trigger objects
   for ( int j = 0; j < tos->size() ; ++j )
   {
      TriggerObject * to = &tos->at(j);
      if ( tonmin_[idx] < 1 ) // all original objects are selected, no emulation
      {
         seltos.push_back(to);
      }
      else  // emulation
      {
         seltos.push_back(to);
      }
   }
   return seltos;

}

std::vector<L1TJet *> SelectTriggerObjects(const std::shared_ptr< Collection<L1TJet> > & tos)
{
   std::vector<L1TJet *> seltos;
   
   int n = tos->size();
   if ( tos->size() > 12 ) n = 12; // L1 decisions are made in the first 12 objects
   
   // select trigger objects
   for ( int j = 0; j < n ; ++j )
   {
      L1TJet * to = &tos->at(j);
      if ( l1tjetsnmin_ < 1 ) // all original objects are selected, no emulation
      {
         seltos.push_back(to);
      }
      else  // emulation
      {
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


bool IsGoodPrimaryVertex(const Vertex & pv)
{
   bool good = false;
   
   if ( ! pv.fake() && pv.ndof() > 4 && fabs(pv.z()) <= 24 && pv.rho() <=2 ) good = true;
   else good = false;
   
   return good;
}

void CreateHistograms(std::map<std::string, TH1F*> & h1, const std::string & label, const int & nmin )
{
   h1[Form("n_%s",label.c_str())] = new TH1F(Form("n_%s",label.c_str())   , "" ,   20,  0   ,  20  );
   for ( int i = 0; i < nmin ; ++i )
   {
      h1[Form("pt%d_%s",i,label.c_str())]  = new TH1F(Form("pt%d_%s",i,label.c_str())  , "" , 500,  0   ,  1000  );
      h1[Form("eta%d_%s",i,label.c_str())] = new TH1F(Form("eta%d_%s",i,label.c_str()) , "" , 100, -2.5 ,   2.5);
      h1[Form("phi%d_%s",i,label.c_str())] = new TH1F(Form("phi%d_%s",i,label.c_str()) , "" , 128, -3.2 ,   3.2);
   }
}

void FillHistograms(std::map<std::string, TH1F*> & h1, const std::string & label, const std::vector<TriggerObject *> & tos, const int & nmin )
{
   h1[Form("n_%s",label.c_str())] -> Fill(tos.size());
   if ( tos.size() < 1 ) return;
   for ( int i = 0; i < nmin ; ++i )
   {
      h1[Form("pt%d_%s",i,label.c_str())]  -> Fill(tos[i]->pt());
      h1[Form("eta%d_%s",i,label.c_str())] -> Fill(tos[i]->eta());
      h1[Form("phi%d_%s",i,label.c_str())] -> Fill(tos[i]->phi());
   }
}

void FillHistograms(std::map<std::string, TH1F*> & h1, const std::string & label, const std::vector<L1TJet *> & tos, const int & nmin )
{
   h1[Form("n_%s",label.c_str())] -> Fill(tos.size());
   if ( tos.size() < 1 ) return;
   for ( int i = 0; i < nmin ; ++i )
   {
      h1[Form("pt%d_%s",i,label.c_str())]  -> Fill(tos[i]->pt());
      h1[Form("eta%d_%s",i,label.c_str())] -> Fill(tos[i]->eta());
      h1[Form("phi%d_%s",i,label.c_str())] -> Fill(tos[i]->phi());
   }
   
}

