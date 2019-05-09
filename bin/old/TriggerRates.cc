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
bool TriggerAccept(Analysis & analysis);

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   // Trigger path info
   analysis.triggerResults(Form("%s/%s",treePath_.c_str(),triggerCol_.c_str()));
   // trigger objects
   for ( auto & obj : triggerObjects_ ) analysis.addTree<TriggerObject> (obj,Form("%s/%s/%s",treePath_.c_str(),triggerObjDir_.c_str(),obj.c_str()));
   
   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
//   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   // counting run by run
   std::map<int,int> totalEvents;
   std::map<int,int> triggerEvents;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%10000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      int run = analysis.run();
      if ( totalEvents.find(run) == totalEvents.end() )  totalEvents[run] = 0;
      if ( triggerEvents.find(run) == totalEvents.end() )  triggerEvents[run] = 0;
      
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      ++totalEvents[run];
      
// Offline selection
// =================
      
      
      
// Trigger analysis, reference trigger (denominator)
// =================================================
      bool trigAccept = TriggerAccept(analysis);
      if ( ! trigAccept ) continue;
      ++triggerEvents[run];
      
   } // end of event loop
   
//    for (auto & ih1 : h1)
//    {
//       ih1.second -> Write();
//    }
//    
//    hout.Close();
   
   int tot = 0;
   int trg = 0;
   for ( auto x: totalEvents )
   {
      std::cout << x.first << "   " << x.second << "   " << triggerEvents[x.first] << std::endl;
      tot += x.second;
      trg += triggerEvents[x.first];
   }
   std::cout << tot << "   " << trg << std::endl;
   
//   std::cout << "TriggerRates: program finished" << std::endl;
   
} //end main


// ================================

bool TriggerAccept(Analysis & analysis)
{
   
// Trigger
// =================================================      
   bool hltAccept = analysis.triggerResult(hltPath_);
   bool l1Accept = true;
   if ( l1Seed_ != "" )  l1Accept = analysis.triggerResult(l1Seed_);
   if ( !hltAccept || !l1Accept ) return false;
   
   // selection of trigger objects
   if ( !trigemul_ ) return true;
   
   for ( size_t i = 0; i < triggerObjects_.size() ; ++i )
   {
      auto to = analysis.collection<TriggerObject>(triggerObjects_[i]);
      auto selto = SelectTriggerObjects(to,i);
      if ( (int)selto.size() < tonmin_[i] ) return false;
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
      else  // emulation for symmetric cuts for multiple objects
      {
         if ( to -> pt() >= toptmin_[idx][0] && fabs( to -> eta() ) <= toetamax_[idx][0] ) 
         {
            seltos.push_back(to);
         }
      }
   }
   return seltos;

}
