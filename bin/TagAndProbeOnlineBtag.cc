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
#include "TF1.h"

#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

using namespace std;
using namespace analysis;
using namespace analysis::trigger;

// =============================================================================================
int main(int argc, char **argv)
{
   TH1::SetDefaultSumw2(); // proper treatment of errors when scaling histograms

   OnlineBtagAnalyser onlinebtag(argc, argv);

   auto config = onlinebtag.config();
   auto workflow = config->workflow();

   // Jets in the analysis
   int tag   = config->tagJet();
   int probe = config->probeJet();

   onlinebtag.jetHistograms(config->nJetsMin(),"tag");
   onlinebtag.jetHistograms(config->nJetsMin(),"probe");
   onlinebtag.jetHistograms(config->nJetsMin(),"probe_match");

   onlinebtag.tagAndProbeTree();

   if ( workflow < 4 )
   {
      for (int i = 0; i < onlinebtag.nEvents(); ++i)
      {
         // ========   GENERAL SELECTION    ==========
         if (!onlinebtag.event(i))                             continue;
         if (!onlinebtag.selectionHLT())                       continue;
         if (!onlinebtag.selectionL1 ())                       continue;
         if (!onlinebtag.selectionPrimaryVertex())             continue;
         if (!onlinebtag.selectionJetId())                     continue;
         if (!onlinebtag.selectionJetPileupId())               continue;
         if (!onlinebtag.selectionNJets())                     continue;
         onlinebtag.actionApplyBjetRegression();                           
         onlinebtag.actionApplyJER();                                      
         if (!onlinebtag.selectionJetDphi(tag, probe))         continue;
         if ( workflow == 2 )
         {
            if (!onlinebtag.dijetTagAndProbe())                continue;
         }

         // ========   TAG SELECTION    ==========
         if (!onlinebtag.selectionJet(tag))                    continue;   
         if (!onlinebtag.onlineJetMatching(tag))               continue;
         if (!onlinebtag.selectionBJet(tag))                   continue;
         onlinebtag.actionApplyBtagSF(tag);
         if (!onlinebtag.onlineBJetMatching(tag))              continue;
         // additional tag selection
         if (!onlinebtag.selectionJetQGlikelihood(tag))        continue;
         if (!onlinebtag.selectionBJetProbBB(tag))             continue;
         if (!onlinebtag.selectionBJetProbC(tag))              continue;
         onlinebtag.fillJetHistograms("tag");

         // ========   PROBE SELECTION    ==========
         if (!onlinebtag.selectionJet(probe))                  continue;
         if (!onlinebtag.onlineJetMatching(probe))             continue;
         if (!onlinebtag.selectionBJet(probe))                 continue;
         onlinebtag.actionApplyBtagSF(probe);
         onlinebtag.fillJetHistograms("probe");
         if ( workflow == 3 )
         {
            onlinebtag.probeMatch(onlinebtag.onlineBJetMatching(probe));
            onlinebtag.fillTagAndProbeTree();
         }
         // ========   PROBE MATCHING    ==========
         if ( workflow != 3)
         {
            if (!onlinebtag.onlineBJetMatching(probe))            continue;
         }
         onlinebtag.fillJetHistograms("probe_match");

      } // end loop
   }

   if ( workflow == 4 )
   {
      for (int i = 0; i < onlinebtag.nEvents(); ++i)
      {
         // ========   GENERAL SELECTION    ==========
         if (!onlinebtag.event(i))                             continue;
         if (!onlinebtag.selectionHLT())                       continue;
         if (!onlinebtag.selectionL1 ())                       continue;
         if (!onlinebtag.selectionPrimaryVertex())             continue;
         if (!onlinebtag.selectionJetId())                     continue;
         if (!onlinebtag.selectionJetPileupId())               continue;
         if (!onlinebtag.selectionNJets())                     continue;
         onlinebtag.actionApplyBjetRegression();                           
         onlinebtag.actionApplyJER();                                  

         // ========   TAG SELECTION    ==========
         if (!onlinebtag.selectionJet(tag))                    continue;
         if (!onlinebtag.onlineJetMatching(tag))               continue;  
         // ========   PROBE SELECTION    ==========
         if (!onlinebtag.selectionJet(probe))                  continue;
         if (!onlinebtag.onlineJetMatching(probe))             continue;

         // =======
         if (!onlinebtag.selectionJetDphi(tag, probe))         continue;

         onlinebtag.tagMatch(onlinebtag.onlineBJetMatching(tag));
         onlinebtag.probeMatch(onlinebtag.onlineBJetMatching(probe));

         // Ntuple
         onlinebtag.fillTagAndProbeTree();


      } // end loop
   }

}