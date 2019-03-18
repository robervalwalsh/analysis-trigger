#ifndef Analysis_Trigger_OnlineBtagAnalyser_h
#define Analysis_Trigger_OnlineBtagAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/Trigger
// Class:      Analysis
// 
/**\class Analysis OnlineBtagAnalyser.cc Analysis/Trigger/src/OnlineBtagAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh
//         Created:  Mon, 05 Mar 2019 18:30:00 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
// 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h" 
// user include files

#include "Analysis/Tools/interface/Analyser.h"

//
// class declaration
//

namespace analysis {
   namespace trigger {

      class OnlineBtagAnalyser : public analysis::tools::Analyser {
         
         public:
            OnlineBtagAnalyser();
            OnlineBtagAnalyser(int argc, char ** argv);
           ~OnlineBtagAnalyser();
           
            // ----------member data ---------------------------
         protected:
            bool do_tree_;
            bool use_ai_;
            std::shared_ptr<TTree> tagprobe_tree_;
            int tag_flavour_;
            int probe_flavour_;
            float ptimbalance_;
            float dphi_;
            float deta_;
            float dr_;
            float tag_qglikelihood_;
            float tag_prob_b_;
            float tag_prob_bb_;
            float tag_prob_lepb_;
            float tag_prob_c_;
            float tag_prob_g_;
            float tag_prob_light_;
            float tag_pt_;
            float tag_eta_;
            float probe_pt_;
            float probe_eta_;
            float tpweight_;
            TMVA::Reader *reader_ai_; 
            std::shared_ptr<TMVA::Tools> tools_ai_;
            TString dir_ai_;
            TString prefix_ai_;
            TString methodName_ai_;
            TString weightfile_ai_;
            std::map<std::string,float> varsf_ai_;
            std::map<std::string,int>   varsi_ai_;
            Int_t spec1_, spec2_;
                           
         private:
            void variables_ai(const int &, const int &);
               
         public:
           virtual bool muonJet(const bool & swap = true);
           void fillTagAndProbeTree();
           void tagAndProbeTree();
           bool selectionAI();

      };
   }
}

#endif  // Analysis_Trigger_OnlineBtagAnalyser_h
