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
   
   bool useJSON = (!isMC_ && json_ != "no_json.txt");

   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Trigger path info
   analysis.triggerResults("MssmHbb/Events/TriggerResults");

   // JSON for data   
   if ( useJSON ) analysis.processJsonFile(json_);
   
   int triggerCounter = 0;
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      if (! isMC_ )
      {
         if ( useJSON )
            if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      bool triggerFired = true;
      if ( hltPathsLogic_ == "OR" ) triggerFired = false;
      for ( auto & path : hltPaths_ )
      {
         if ( hltPathsLogic_ == "OR" )
            triggerFired = (triggerFired || analysis.triggerResult(path));
         else
            triggerFired = (triggerFired && analysis.triggerResult(path));
      }
      
      bool triggerFired2 = true;
      if ( hltPathsLogic2_ == "OR" ) triggerFired2 = false;
      for ( auto & path : hltPaths2_ )
      {
         if ( hltPathsLogic2_ == "OR" )
            triggerFired2 = (triggerFired2 || analysis.triggerResult(path));
         else
            triggerFired2 = (triggerFired2 && analysis.triggerResult(path));
      }
      
      if ( triggerFired && !triggerFired2 )
         ++triggerCounter;
      
      
   }
   
// PRINT OUTS   

   std::cout << "Number of pure counts in " << nlumis_ << " = " << triggerCounter << std::endl;
   std::cout << "Pure rate = " << triggerCounter/(nlumis_*23.4) << " Hz" << std::endl;
   
} //end main

