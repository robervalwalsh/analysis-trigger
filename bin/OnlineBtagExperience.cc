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

//#include "Analysis/Tools/interface/Analyser.h"
#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

using namespace std;
using namespace analysis;
//using namespace analysis::tools;
using namespace analysis::trigger;

using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;
     

// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   OnlineBtagAnalyser onlinebtag(argc,argv);
   
   onlinebtag.jetHistograms(2,"tag");
   onlinebtag.jetHistograms(2,"probe");
   onlinebtag.jetHistograms(2,"probe_match");
   
   onlinebtag.tagAndProbeTree();
   
   // Analysis of events
   std::cout << "The sample size is " << onlinebtag.analysis()->size() << " events" << std::endl;
   std::cout << "---------------------------" << std::endl;
   
   std::cout << "Workflow index = " << onlinebtag.config()->workflow() << std::endl;
   std::cout << "--------------------" << std::endl;
   
   int seed = onlinebtag.seed();
   
// 
   for ( int i = 0 ; i < onlinebtag.nEvents() ; ++i )
   {
      bool goodEvent = onlinebtag.event(i);

      if ( ! goodEvent ) continue;

      if ( onlinebtag.config()->workflow() == 1 )  // ========== DATA and MC with data-like sequence ========
      {
      // trigger selection
         if ( ! onlinebtag.selectionHLT()                 )   continue;
         if ( ! onlinebtag.selectionL1 ()                 )   continue;  // to be used in case of "OR" of seeds
      
            
      // jet identification selection
         if ( ! onlinebtag.selectionJetId()               )   continue;
         if ( ! onlinebtag.selectionJetPileupId()         )   continue;
         if ( ! onlinebtag.selectionNJets()               )   continue;
         
         // jet energy resolution  
         onlinebtag.actionApplyJER();
         
         // pileup weight
         onlinebtag.actionApplyPileupWeight();
      
      // delta Phi
         if ( ! onlinebtag.selectionJetDphi(1,2)          )   continue;
      // pT imbalance
         if ( ! onlinebtag.selectionJetPtImbalance(1,2)   )   continue;
         
      // tag jet selection
         if ( ! onlinebtag.selectionJet(2)                )   continue;
         if ( ! onlinebtag.onlineJetMatching(2)           )   continue;
         if ( ! onlinebtag.selectionBJet(2)               )   continue;
         // btag SF
         onlinebtag.actionApplyBtagSF(2);
         if ( ! onlinebtag.onlineBJetMatching(2)          )   continue;
      // additional tag selection
         if ( ! onlinebtag.selectionJetQGlikelihood(2)    )   continue;
         if ( ! onlinebtag.selectionBJetProbBB(2)         )   continue;
         if ( ! onlinebtag.selectionBJetProbC(2)          )   continue;
         
         onlinebtag.fillJetHistograms("tag");
         
      // probe jet selection
      // qg likelihood
         if ( ! onlinebtag.selectionJet(1)                )   continue;
         if ( ! onlinebtag.selectionJetQGlikelihood(1)    )   continue;
         if ( ! onlinebtag.onlineJetMatching(1)           )   continue;
         if ( ! onlinebtag.selectionBJet(1)               )   continue;
         // btag SF
         onlinebtag.actionApplyBtagSF(1);
         
         
      // If there is any AI selection, otherwise always true
         if ( ! onlinebtag.selectionAI()                  )   continue;
      // fill tree for AI
         onlinebtag.fillTagAndProbeTree();

      // fill tag histograms after selection
         onlinebtag.fillJetHistograms("probe");
         
      // PROBE jet match to online
         if ( ! onlinebtag.onlineBJetMatching(1)          )   continue;
         
      // fill probe histograms after selection
         onlinebtag.fillJetHistograms("probe_match");
      }
      
      
//       if ( onlinebtag.config()->workflow() == 2 )  // ========== DATA and MC with data-like sequence - SEMILEPTONIC ========
//       {
//       
//       // trigger selection
//          if ( ! onlinebtag.selectionHLT()                 )   continue;
//          if ( ! onlinebtag.selectionL1 ()                 )   continue;  // to be used in case of "OR" of seeds
//       
//       // muon identification selection
//          if ( ! onlinebtag.selectionMuonId()                 )   continue;
//          if ( ! onlinebtag.selectionNMuons()                 )   continue;
//          if ( ! onlinebtag.selectionMuons()                  )   continue;
//          if ( ! onlinebtag.onlineMuonMatching()              )   continue; // FIXME there: should change the selectedMuons ????
//             
//       // jet identification selection
//          if ( ! onlinebtag.selectionJetId()               )   continue;
//          if ( ! onlinebtag.selectionJetPileupId()         )   continue;
//          if ( ! onlinebtag.selectionNJets()               )   continue;
//          
//       // delta Phi
//          if ( ! onlinebtag.selectionJetDphi(1,2)          )   continue;
//       // pT imbalance
//          if ( ! onlinebtag.selectionJetPtImbalance(1,2)   )   continue;
//          
//       // Muon jet association
//          if ( ! onlinebtag.muonJet()                      )   continue;
//          
//       // tag jet selection
//          if ( ! onlinebtag.selectionJet(2)                )   continue;
//          if ( ! onlinebtag.selectionJetQGlikelihood(2)    )   continue;
//          if ( ! onlinebtag.onlineJetMatching(2)           )   continue;
//          if ( ! onlinebtag.selectionBJet(2)               )   continue;
//          if ( ! onlinebtag.onlineBJetMatching(2)          )   continue;
//          onlinebtag.fillJetHistograms("tag");
//          
//       // probe jet selection
//       // qg likelihood
//          if ( ! onlinebtag.selectionJet(1)                )   continue;
//          if ( ! onlinebtag.selectionJetQGlikelihood(1)    )   continue;
//          if ( ! onlinebtag.onlineJetMatching(1)           )   continue;
//          if ( ! onlinebtag.selectionBJet(1)               )   continue;
//          
//       // fill tag histograms after selection
//          onlinebtag.fillJetHistograms("probe");
//          
//       // PROBE jet match to online
//          if ( ! onlinebtag.onlineBJetMatching(1)          )   continue;
//          
//       // fill probe histograms after selection
//          onlinebtag.fillJetHistograms("probe_match");
//       }
      
   }
   
} //end main

