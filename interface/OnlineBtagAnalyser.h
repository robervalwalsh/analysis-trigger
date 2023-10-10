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
            /// General parameters
            int   process_id_;
            /// Correlations
            float ptimbalance_;
            float dphi_;
            float deta_;
            float dr_;            
            /// Probe variables
            int probe_jet_;
            int probe_flavour_;
            float probe_qglikelihood_;
            float probe_prob_b_;
            float probe_prob_bb_;
            float probe_prob_lepb_;
            float probe_prob_c_;
            float probe_prob_g_;
            float probe_prob_light_;
            float probe_pt_;
            float probe_eta_;
            float probe_phi_;
            int   probe_match_;
            float probe_btag_sfl_;
            float probe_btag_sfm_;
            float probe_btag_sft_;
            /// tag variables
            int tag_jet_;
            int tag_flavour_;
            float tag_qglikelihood_;
            float tag_prob_b_;
            float tag_prob_bb_;
            float tag_prob_lepb_;
            float tag_prob_c_;
            float tag_prob_g_;
            float tag_prob_light_;
            float tag_pt_;
            float tag_eta_;
            float tag_phi_;
            int   tag_match_;
            float tag_btag_sfl_;
            float tag_btag_sfm_;
            float tag_btag_sft_;
            /// Jet3 varibles
            float jet3_pt_;
            float jet3_eta_;
            float jet3_phi_;
            /// event weight 
            float tpweight_;
            /// primary vertices
            int num_primary_vertices_;
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
            void variables_ai();
               
         public:
           virtual bool muonJet(const bool & swap = true);
           void fillTagAndProbeTree();
           void tagAndProbeTree();
           bool selectionAI();
           /// veto events with a 3rd jet
           bool dijetTagAndProbe();
           /// probe match
           void probeMatch(const bool &);
           /// tag match
           void tagMatch(const bool &);

      };
   }
}

#endif  // Analysis_Trigger_OnlineBtagAnalyser_h
