/**\class Trigger OnlineBtagAnalyser.cc Analysis/Trigger/src/OnlineBtagAnalyser.cc

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
#include <iostream>
// user include files
#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::trigger;
using namespace TMVA;


//
// constructors and destructor
//
OnlineBtagAnalyser::OnlineBtagAnalyser()
{
}

OnlineBtagAnalyser::OnlineBtagAnalyser(int argc, char ** argv) : BaseAnalyser(argc,argv), Analyser(argc,argv)
{
   do_tree_ = config_->doTree();
   use_ai_ = false;
   reader_ai_ = nullptr;
   tools_ai_  = nullptr;
   /// General parameters
   process_id_ = -1;
   /// Correlations
   ptimbalance_ = -1;
   dphi_ = -1;
   deta_ = -1;
   dr_ = -1;
   /// Probe variables
   probe_jet_ = config_->probeJet()-1;
   probe_flavour_ = -1;
   probe_qglikelihood_ = -1;
   probe_prob_b_ = -1;
   probe_prob_bb_ = -1;
   probe_prob_lepb_ = -1;
   probe_prob_c_ = -1;
   probe_prob_g_ = -1;
   probe_prob_light_ = -1;
   probe_pt_ = -1;
   probe_eta_ = -10;
   probe_phi_ = -100;
   probe_match_ = 0;
   probe_btag_sfl_ = 1;
   probe_btag_sfm_ = 1;
   probe_btag_sft_ = 1;

   /// tag variables
   tag_jet_ = config_->tagJet()-1;
   tag_flavour_ = -1;
   tag_qglikelihood_ = -1;
   tag_prob_b_ = -1;
   tag_prob_bb_ = -1;
   tag_prob_lepb_ = -1;
   tag_prob_c_ = -1;
   tag_prob_g_ = -1;
   tag_prob_light_ = -1;
   tag_pt_ = -1;
   tag_eta_ = -10;
   tag_phi_ = -100;
   tag_match_ = 0;
   tag_btag_sfl_ = 1;
   tag_btag_sfm_ = 1;
   tag_btag_sft_ = 1;


   /// Jet3 variables
   jet3_pt_ = -1;
   jet3_eta_ = -10;
   jet3_phi_ = -100;
   /// event weight
   tpweight_ = 0;
   
   if ( config_->methodAI() != "" || config_->directoryAI() != "" || config_->variablesAI("F").size() != 0 )
   {
      use_ai_ = true;
      tools_ai_ = std::shared_ptr<TMVA::Tools>(&(TMVA::Tools::Instance()));
      // Create the Reader object
      reader_ai_     = new TMVA::Reader( "!Color:!Silent" );
      dir_ai_        = TString(config_ -> directoryAI()+"/");
      prefix_ai_     = "TMVAClassification";
      methodName_ai_ = TString(config_->methodAI()) + TString(" method");
      weightfile_ai_ = dir_ai_ + prefix_ai_ + TString("_") + TString(config_->methodAI()) + TString(".weights.xml");
      
      reader_ai_->AddSpectator( "probe_flavour", &spec1_ );
      reader_ai_->AddSpectator( "tag_flavour",   &spec2_ );
      
      varsf_ai_["ptimbalance"]      = ptimbalance_;
      varsf_ai_["dphi"]             = dphi_;
      varsf_ai_["deta"]             = deta_;
      varsf_ai_["dr"]               = dr_;
      varsf_ai_["tag_qglikelihood"] = tag_qglikelihood_;
      varsf_ai_["tag_prob_b"]       = tag_prob_b_;
      varsf_ai_["tag_prob_bb"]      = tag_prob_bb_;
      varsf_ai_["tag_prob_lepb"]    = tag_prob_lepb_;
      varsf_ai_["tag_prob_c"]       = tag_prob_c_;
      varsf_ai_["tag_prob_g"]       = tag_prob_g_;
      varsf_ai_["tag_prob_light"]   = tag_prob_light_;
      varsf_ai_["tag_pt"]           = tag_pt_;
      varsf_ai_["abs(tag_eta)"]     = -fabs(tag_eta_);
      varsf_ai_["probe_pt"]         = probe_pt_;
      varsf_ai_["abs(probe_eta)"]   = -fabs(probe_eta_);
      for ( auto & var : config_ -> variablesAI("F") )
      {
         TString v(var);
         reader_ai_->AddVariable( v, &varsf_ai_[var] );
      }
   
   // Book the MVA methods

   reader_ai_->BookMVA( methodName_ai_, weightfile_ai_ );
   

   
   }

}

OnlineBtagAnalyser::~OnlineBtagAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if ( do_tree_ ) tagprobe_tree_ -> Write();
}


//
// member functions
//
// ------------ method called for each event  ------------


bool OnlineBtagAnalyser::muonJet(const bool & swap)
{
   // jet with ranking 1 is the muon jet, swap with jet 2 in case it has the leading muon 
   int r1 = 1;
   int r2 = 2;
   int j1 = r1-1;
   int j2 = r2-1;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
   {
      if ( swap ) h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Online Btag T&P Semileptonic: Jet-muon association -> Muon-Jet index 1 (probe jet)");
      else        h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Online Btag T&P Semileptonic: Jet-muon association");
   }
   
   auto jet1 = selectedJets_[j1];
   jet1 -> addMuon(selectedMuons_);
   auto jet2 = selectedJets_[j2];
   jet2 -> addMuon(selectedMuons_);
   
   if ( ! (jet1 -> muon() || jet2 -> muon()) ) return false;
   if ( !  jet1 -> muon() && swap ) this->jetSwap(r1,r2);
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);
   return true;
   
}

void OnlineBtagAnalyser::tagAndProbeTree()
{
   if ( ! do_tree_ ) return;
   this->output()->cd();
   tagprobe_tree_ = std::make_shared<TTree>("TagAndProbe","TTree with variables for the online btag tag and probe");
   tagprobe_tree_ -> Branch("process_id",&process_id_,"process_id/I");
   /// Correlations
   tagprobe_tree_ -> Branch("ptimbalance",&ptimbalance_,"ptimbalance/F");
   tagprobe_tree_ -> Branch("dphi",&dphi_,"dphi/F");
   tagprobe_tree_ -> Branch("deta",&deta_,"deta/F");
   tagprobe_tree_ -> Branch("dr",&dr_,"dr/F");

   /// Probe variables
   tagprobe_tree_ -> Branch("probe_jet",&probe_jet_,"probe_jet/I");
   tagprobe_tree_ -> Branch("probe_flavour",&probe_flavour_,"probe_flavour/I");
   tagprobe_tree_ -> Branch("probe_qglikelihood",&probe_qglikelihood_,"probe_qglikelihood/F");
   tagprobe_tree_ -> Branch("probe_prob_b",&probe_prob_b_,"probe_prob_b/F");
   tagprobe_tree_ -> Branch("probe_prob_bb",&probe_prob_bb_,"probe_prob_bb/F");
   tagprobe_tree_ -> Branch("probe_prob_lepb",&probe_prob_lepb_,"probe_prob_lepb/F");
   tagprobe_tree_ -> Branch("probe_prob_c",&probe_prob_c_,"probe_prob_c_/F");
   tagprobe_tree_ -> Branch("probe_prob_g",&probe_prob_g_,"probe_prob_g/F");
   tagprobe_tree_ -> Branch("probe_prob_light",&probe_prob_light_,"probe_prob_light/F");
   tagprobe_tree_ -> Branch("probe_pt",&probe_pt_,"probe_pt/F");
   tagprobe_tree_ -> Branch("probe_eta",&probe_eta_,"probe_eta/F");
   tagprobe_tree_ -> Branch("probe_phi",&probe_phi_,"probe_phi/F");
   tagprobe_tree_ -> Branch("probe_match",&probe_match_,"probe_match/I");
   tagprobe_tree_ -> Branch("probe_btag_sfl",&probe_btag_sfl_,"probe_btag_sfl/F");
   tagprobe_tree_ -> Branch("probe_btag_sfm",&probe_btag_sfm_,"probe_btag_sfm/F");
   tagprobe_tree_ -> Branch("probe_btag_sft",&probe_btag_sft_,"probe_btag_sft/F");

   /// Tag variables
   tagprobe_tree_ -> Branch("tag_jet",&tag_jet_,"tag_jet/I");
   tagprobe_tree_ -> Branch("tag_flavour",&tag_flavour_,"tag_flavour/I");
   tagprobe_tree_ -> Branch("tag_qglikelihood",&tag_qglikelihood_,"tag_qglikelihood/F");
   tagprobe_tree_ -> Branch("tag_prob_b",&tag_prob_b_,"tag_prob_b/F");
   tagprobe_tree_ -> Branch("tag_prob_bb",&tag_prob_bb_,"tag_prob_bb/F");
   tagprobe_tree_ -> Branch("tag_prob_lepb",&tag_prob_lepb_,"tag_prob_lepb/F");
   tagprobe_tree_ -> Branch("tag_prob_c",&tag_prob_c_,"tag_prob_c_/F");
   tagprobe_tree_ -> Branch("tag_prob_g",&tag_prob_g_,"tag_prob_g/F");
   tagprobe_tree_ -> Branch("tag_prob_light",&tag_prob_light_,"tag_prob_light/F");
   tagprobe_tree_ -> Branch("tag_pt",&tag_pt_,"tag_pt/F");
   tagprobe_tree_ -> Branch("tag_eta",&tag_eta_,"tag_eta/F");
   tagprobe_tree_ -> Branch("tag_phi",&tag_phi_,"tag_phi/F");
   tagprobe_tree_ -> Branch("tag_match",&tag_match_,"tag_match/I");
   tagprobe_tree_ -> Branch("tag_btag_sfl",&tag_btag_sfl_,"tag_btag_sfl/F");
   tagprobe_tree_ -> Branch("tag_btag_sfm",&tag_btag_sfm_,"tag_btag_sfm/F");
   tagprobe_tree_ -> Branch("tag_btag_sft",&tag_btag_sft_,"tag_btag_sft/F");

   /// Jet3 variables
   tagprobe_tree_ -> Branch("jet3_pt",&jet3_pt_,"jet3_pt/F");
   tagprobe_tree_ -> Branch("jet3_eta",&jet3_eta_,"jet3_eta/F");
   tagprobe_tree_ -> Branch("jet3_phi",&jet3_phi_,"jet3_phi/F");
   // event weight 
   tagprobe_tree_ -> Branch("weight",&tpweight_,"weight/F");
   
}

void OnlineBtagAnalyser::fillTagAndProbeTree()
{
   if ( ! do_tree_ ) return;
   std::string algo = config_->btagAlgorithm();
   if ( ! do_tree_ || algo != "deepflavour" ) return;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Fill TagAndProbe tree");

   // variables for the tmva tree
   variables_ai();   
   
   tpweight_         = weight_;
   
   tagprobe_tree_ -> Fill();
   h1_["cutflow"] -> Fill(cutflow_,weight_);

}

void OnlineBtagAnalyser::variables_ai()
{
   process_id_       = config_->processId();
   ptimbalance_      = fabs(selectedJets_[probe_jet_]->pt() - selectedJets_[tag_jet_]->pt())/selectedJets_[probe_jet_]->pt();
   dphi_             = fabs(selectedJets_[probe_jet_]->deltaPhi(*selectedJets_[tag_jet_]));
   deta_             = fabs(selectedJets_[probe_jet_]->eta() - selectedJets_[tag_jet_]->eta());
   dr_               = fabs(selectedJets_[probe_jet_]->deltaR(*selectedJets_[tag_jet_]));

   probe_pt_         = selectedJets_[probe_jet_]->pt();
   probe_eta_        = selectedJets_[probe_jet_]->eta();
   probe_phi_        = selectedJets_[probe_jet_]->phi();
   probe_qglikelihood_ = selectedJets_[probe_jet_]->qgLikelihood();
   probe_prob_b_       = selectedJets_[probe_jet_]->btag("btag_dfb");
   probe_prob_bb_      = selectedJets_[probe_jet_]->btag("btag_dfbb");
   probe_prob_lepb_    = selectedJets_[probe_jet_]->btag("btag_dflepb");
   probe_prob_c_       = selectedJets_[probe_jet_]->btag("btag_dfc");
   probe_prob_g_       = selectedJets_[probe_jet_]->btag("btag_dfg");
   probe_prob_light_   = selectedJets_[probe_jet_]->btag("btag_dflight");


   tag_qglikelihood_ = selectedJets_[tag_jet_]->qgLikelihood();
   tag_prob_b_       = selectedJets_[tag_jet_]->btag("btag_dfb");
   tag_prob_bb_      = selectedJets_[tag_jet_]->btag("btag_dfbb");
   tag_prob_lepb_    = selectedJets_[tag_jet_]->btag("btag_dflepb");
   tag_prob_c_       = selectedJets_[tag_jet_]->btag("btag_dfc");
   tag_prob_g_       = selectedJets_[tag_jet_]->btag("btag_dfg");
   tag_prob_light_   = selectedJets_[tag_jet_]->btag("btag_dflight");
   tag_pt_           = selectedJets_[tag_jet_]->pt();
   tag_eta_          = selectedJets_[tag_jet_]->eta();
   tag_phi_          = selectedJets_[tag_jet_]->phi();


   /// Jet3 variables
   if ( selectedJets_.size() > 2 )
   {
      jet3_pt_       = selectedJets_[2]->pt();
      jet3_eta_      = selectedJets_[2]->eta();
      jet3_phi_      = selectedJets_[2]->phi();
   }

   /// Flavours of tag and probe 
   if ( config_->isMC() && config_ -> histogramJetsPerFlavour() )
   {
      if ( ! config_ -> useJetsExtendedFlavour() )
      {
         probe_flavour_ = abs(selectedJets_[probe_jet_]->flavour());
         tag_flavour_   = abs(selectedJets_[tag_jet_]->flavour());
      }
      else
      {
         std::string probe_flv = selectedJets_[probe_jet_]->extendedFlavour();
         std::string tag_flv   = selectedJets_[tag_jet_]->extendedFlavour();
         if      ( probe_flv == "b"  ) probe_flavour_ = 5;
         else if ( probe_flv == "c"  ) probe_flavour_ = 4;
         else if ( probe_flv == "bb" ) probe_flavour_ = 55;
         else if ( probe_flv == "cc" ) probe_flavour_ = 44;
         else                          probe_flavour_ = 0;

         if      ( tag_flv == "b"    ) tag_flavour_ = 5;
         else if ( tag_flv == "c"    ) tag_flavour_ = 4;
         else if ( tag_flv == "bb"   ) tag_flavour_ = 55;
         else if ( tag_flv == "cc"   ) tag_flavour_ = 44;
         else                          tag_flavour_ = 0; 

      }
   }
   /// btag scale factors
   if ( config_->isMC() )
   {
      probe_btag_sfl_ = this->btagSF(probe_jet_+1, "loose").nominal;
      probe_btag_sfm_ = this->btagSF(probe_jet_+1, "medium").nominal;
      probe_btag_sft_ = this->btagSF(probe_jet_+1, "tight").nominal;
      tag_btag_sfl_ = this->btagSF(tag_jet_+1, "loose").nominal;
      tag_btag_sfm_ = this->btagSF(tag_jet_+1, "medium").nominal;
      tag_btag_sft_ = this->btagSF(tag_jet_+1, "tight").nominal;
   }
   
   varsf_ai_["ptimbalance"]      = ptimbalance_;
   varsf_ai_["dphi"]             = dphi_;
   varsf_ai_["deta"]             = deta_;
   varsf_ai_["dr"]               = dr_;
   varsf_ai_["tag_qglikelihood"] = tag_qglikelihood_;
   varsf_ai_["tag_prob_b"]       = tag_prob_b_;
   varsf_ai_["tag_prob_bb"]      = tag_prob_bb_;
   varsf_ai_["tag_prob_lepb"]    = tag_prob_lepb_;
   varsf_ai_["tag_prob_c"]       = tag_prob_c_;
   varsf_ai_["tag_prob_g"]       = tag_prob_g_;
   varsf_ai_["tag_prob_light"]   = tag_prob_light_;
   varsf_ai_["tag_pt"]           = tag_pt_;
   varsf_ai_["abs(tag_eta)"]     = fabs(tag_eta_);
   varsf_ai_["probe_pt"]         = probe_pt_;
   varsf_ai_["abs(probe_eta)"]   = fabs(probe_eta_);
   
}

bool OnlineBtagAnalyser::selectionAI()
{
   if ( ! reader_ai_ ) return true;
   if ( config_->discriminatorMaxAI() < -1000. && config_->discriminatorMinAI() < -1000. && config_->efficiencyMinAI() < 0. ) return true;
   
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
   {
      if ( TString(config_->methodAI()).BeginsWith("Cuts") )
      {
         if ( config_->efficiencyMinAI() > 0 )
            h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,Form("AI %s: efficiency > %5.3f",methodName_ai_.Data(),config_->efficiencyMinAI()));
      }
      else
      {
         if ( config_->discriminatorMaxAI() > -1000 && config_->discriminatorMinAI() < -1000 )
            h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,Form("AI %s: discriminator < %5.3f",methodName_ai_.Data(),config_->discriminatorMaxAI()));
         if ( config_->discriminatorMinAI() > -1000 && config_->discriminatorMaxAI() < config_->discriminatorMinAI() )
            h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,Form("AI %s: discriminator > %5.3f",methodName_ai_.Data(),config_->discriminatorMinAI()));
         if ( config_->discriminatorMinAI() > -1000 && config_->discriminatorMaxAI() > -1000 )
            h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,Form("AI %s: %5.3f < discriminator < %5.3f",methodName_ai_.Data(),config_->discriminatorMinAI(),config_->discriminatorMaxAI()));
      }
   }
   variables_ai();
   
   // all possible variables, but if training done with less that will be controlled from the config
   // vector<float> varsf_ai_ | vector<int> varsi_ai_  (see loops below)
   
   if ( TString(config_->methodAI()).BeginsWith("Cuts") )
   {
      bool cuts = reader_ai_->EvaluateMVA( methodName_ai_, config_->efficiencyMinAI() );
      if ( ! cuts ) return false;
   }
   else
   {
      double discr = reader_ai_->EvaluateMVA( methodName_ai_ );
      if ( config_->discriminatorMaxAI() > -1000 && discr > config_->discriminatorMaxAI() ) return false;
      if ( config_->discriminatorMinAI() > -1000 && discr < config_->discriminatorMinAI() ) return false;
   }
      
   h1_["cutflow"] -> Fill(cutflow_,weight_);
   return true;
}


bool OnlineBtagAnalyser::dijetTagAndProbe()
{
   bool is_dijet = true;
   std::string label = "3rd jet veto: if pt3 > 0.3*(pt1+pt2)/2";
   if ( selectedJets_.size() > 2 )
   {
      float mean_pt1_pt2 = (selectedJets_[0]->pt()+selectedJets_[1]->pt())/2;
      is_dijet = (selectedJets_[2]->pt() < 0.3*mean_pt1_pt2);
   }

   cutflow(label, is_dijet);

   return is_dijet;
}

void OnlineBtagAnalyser::probeMatch(const bool & match)
{
   if ( match ) probe_match_ = 1;
   else probe_match_ = 0;
}

void OnlineBtagAnalyser::tagMatch(const bool & match)
{
   if ( match ) tag_match_ = 1;
   else tag_match_ = 0;
}