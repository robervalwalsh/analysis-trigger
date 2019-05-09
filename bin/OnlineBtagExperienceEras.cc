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
#include "TRandom3.h"

//#include "Analysis/Tools/interface/Analyser.h"
#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

using namespace std;
using namespace analysis;
//using namespace analysis::tools;
using namespace analysis::trigger;

using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;
     
TH1F * hera_;
bool ProcessEvent(OnlineBtagAnalyser &, const std::string & era = "");
int Eras(OnlineBtagAnalyser &);

// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   OnlineBtagAnalyser onlinebtag(argc,argv);
   
   int eras = Eras(onlinebtag);
   if ( eras < 0 ) return -1;
   
   
   // Analysis of events
   std::cout << "The sample size is " << onlinebtag.analysis()->size() << " events" << std::endl;
   std::cout << "---------------------------" << std::endl;
   
   std::cout << "Workflow index = " << onlinebtag.config()->workflow() << std::endl;
   std::cout << "--------------------" << std::endl;
   
   int seed = onlinebtag.seed();
   TRandom3 * rnd = new TRandom3(seed);
   
// 
   for ( int i = 0 ; i < onlinebtag.nEvents() ; ++i )
   {
      bool goodEvent = onlinebtag.event(i);

      if ( ! goodEvent ) continue;

      if ( onlinebtag.config()->workflow() == 1 )  // ========== DATA and MC with data-like sequence ========
      {
         if ( hera_ )
         {
            auto x = rnd->Rndm();
            std::string era = (onlinebtag.config()->era())[hera_->FindBin(x)-1];

               if ( ! ProcessEvent(onlinebtag,era) ) continue;
         }
         else
         {
            if ( ! ProcessEvent(onlinebtag) ) continue;
         }
      }
      
   }
   
} //end main

// ===================================================================
bool ProcessEvent(OnlineBtagAnalyser & onlinebtag, const std::string & era)
{
   
   std::string hname;
   
// trigger selection
   if ( ! onlinebtag.selectionHLT()                 )   return false;
   if ( ! onlinebtag.selectionL1 ()                 )   return false;  // to be used in case of "OR" of seeds

      
// jet identification selection
   if ( ! onlinebtag.selectionJetId()               )   return false;
   if ( ! onlinebtag.selectionJetPileupId()         )   return false;
   if ( ! onlinebtag.selectionNJets()               )   return false;
   
      
   // jet energy resolution  
   onlinebtag.actionApplyJER();

// delta Phi
   if ( ! onlinebtag.selectionJetDphi(1,2)          )   return false;
// pT imbalance
   if ( ! onlinebtag.selectionJetPtImbalance(1,2)   )   return false;
   
// tag jet selection
   if ( ! onlinebtag.selectionJet(2)                )   return false;
   if ( ! onlinebtag.onlineJetMatching(2)           )   return false;
   if ( ! onlinebtag.selectionBJet(2)               )   return false;
   // btag SF
   onlinebtag.actionApplyBtagSF(2);
   if ( ! onlinebtag.onlineBJetMatching(2)          )   return false;
// additional tag selection
   if ( ! onlinebtag.selectionJetQGlikelihood(2)    )   return false;
   if ( ! onlinebtag.selectionBJetProbBB(2)         )   return false;
   if ( ! onlinebtag.selectionBJetProbC(2)          )   return false;
   
   hname = "tag";
   if ( era != "" ) hname = hname+"_"+era;

   onlinebtag.fillJetHistograms(hname.c_str());
   
// probe jet selection
// qg likelihood
   if ( ! onlinebtag.selectionJet(1)                )   return false;
   if ( ! onlinebtag.selectionJetQGlikelihood(1)    )   return false;
   if ( ! onlinebtag.onlineJetMatching(1)           )   return false;
   if ( ! onlinebtag.selectionBJet(1)               )   return false;
   // btag SF
   onlinebtag.actionApplyBtagSF(1);
   
   
// If there is any AI selection, otherwise always true
   if ( ! onlinebtag.selectionAI()                  )   return false;
// fill tree for AI
   onlinebtag.fillTagAndProbeTree();

// fill tag histograms after selection
   hname = "probe";
   if ( era != "" ) hname = hname+"_"+era;
   onlinebtag.fillJetHistograms(hname.c_str());
   
// PROBE jet match to online
   if ( ! onlinebtag.onlineBJetMatching(1)          )   return false;
   
// fill probe histograms after selection
   hname = "probe_match";
   if ( era != "" ) hname = hname+"_"+era;
   onlinebtag.fillJetHistograms(hname.c_str());


   return true;
}

// ===================================================================
// ===================================================================
// ===================================================================
int Eras(OnlineBtagAnalyser & onlinebtag)
{
   std::vector<std::string> eras    = onlinebtag.config()->era();
   std::vector<float>       eralumi = onlinebtag.config()->eraLumi();
   
   if ( eras.size() != eralumi.size() )
   {
      std::cout << std::endl;
      std::cout << "======================================================" << std::endl;
      std::cout << "Number of eras different from amount of lumi values!!!" << std::endl;
      std::cout << "Check your configuration file!!!" << std::endl;     
      std::cout << "======================================================" << std::endl;
      std::cout << std::endl;
      
      return -1;
   }
   
   if ( eras.size() == 0 )
   {
      onlinebtag.jetHistograms(2,"tag");
      onlinebtag.jetHistograms(2,"probe");
      onlinebtag.jetHistograms(2,"probe_match");
      hera_ = NULL;
      return 0;
   }
      onlinebtag.tagAndProbeTree();
   
   std::vector<float>       fraclumi;
   std::vector<float>       cumlumi;
   float totallumi = 0;
   for ( auto & el: eralumi ) totallumi += el;
   
   cumlumi.push_back(0);
   for ( auto & era : eras )
   {
      auto i = &era - &eras[0];
      float fl = eralumi[i]/totallumi;
      fraclumi.push_back(fl);
      cumlumi.push_back(cumlumi.back()+fl);
      onlinebtag.jetHistograms(2,Form("tag_%s",era.c_str()));
      onlinebtag.jetHistograms(2,Form("probe_%s",era.c_str()));
      onlinebtag.jetHistograms(2,Form("probe_match_%s",era.c_str()));
   }
   cumlumi.back() = 1.;
   hera_ = new TH1F("hera","",cumlumi.size(),&cumlumi[0]);
   return 1;
   
}
