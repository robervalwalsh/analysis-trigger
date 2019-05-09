// Tag & Probe for online btag efficiency
// tag jet is second leading in pT, probe jet is leading in pT
//
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

#include "Analysis/Tools/interface/Analyser.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;

using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;
     

// float GetBTag(const Jet & jet, const std::string & algo);
// void CreateHistograms(TH1s &, const int & n = 3);
// void CreateHistograms(TH2s &, const int & n = 3);
// template <typename T>
// void WriteHistograms(T & );


// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   Analyser onlinebtag(argc,argv);
   
   onlinebtag.jetHistograms(2,"nomatch_1");
   onlinebtag.jetHistograms(2,"match_1");
   
   // Analysis of events
   std::cout << "The sample size is " << onlinebtag.analysis()->size() << " events" << std::endl;
   
// 
   for ( int i = 0 ; i < onlinebtag.nEvents() ; ++i )
   {
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      bool goodEvent = onlinebtag.event(i);
      if ( ! goodEvent ) continue;
      
   // trigger selection
      if ( ! onlinebtag.selectionTrigger()       )   continue;
         
   // jet identification selection
      if ( ! onlinebtag.selectionJetId()         )   continue;
      if ( ! onlinebtag.selectionJetPileupId()   )   continue;
      if ( ! onlinebtag.selectionNJets()         )   continue;
      
   //  1st and 2nd jet kinematic selection
      if ( ! onlinebtag.selectionJet(1)          )   continue;
      if ( ! onlinebtag.selectionJet(2)          )   continue;
      
   // btag of two leading jets
      if ( ! onlinebtag.selectionBJet(1)         )   continue;
      if ( ! onlinebtag.selectionBJet(2)         )   continue;
      
   // delta phi jet selection
      if ( ! onlinebtag.selectionJetDphi(1,2)    )   continue;
      
   // jets 1,2 matching to online jets
      if ( ! onlinebtag.onlineJetMatching(1)     )   continue;
      if ( ! onlinebtag.onlineJetMatching(2)     )   continue;
      
   // jet 2 (TAG) matching to online btag objects
      if ( ! onlinebtag.onlineBJetMatching(2)    )   continue;
      
      onlinebtag.fillJetHistograms("nomatch_1");
           
   // jet 1 (PROBE) matching to online btag objects
      if ( ! onlinebtag.onlineBJetMatching(1)    )   continue;
      
      onlinebtag.fillJetHistograms("match_1");
      
   }
   
} //end main

