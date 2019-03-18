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
   tag_flavour_ = 0;
   probe_flavour_ = 0;
   ptimbalance_ = -1;
   dphi_ = -1;
   deta_ = -1;
   dr_ = -1;
   tag_qglikelihood_ = -1;
   tag_prob_b_ = -1;
   tag_prob_bb_ = -1;
   tag_prob_lepb_ = -1;
   tag_prob_c_ = -1;
   tag_prob_g_ = -1;
   tag_prob_light_ = -1;
   tag_pt_ = -1;
   tag_eta_ = -10;
   
   
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
   tagprobe_tree_ -> Branch("tag_flavour",&tag_flavour_,"tag_flavour/I");
   tagprobe_tree_ -> Branch("probe_flavour",&probe_flavour_,"probe_flavour/I");
   tagprobe_tree_ -> Branch("ptimbalance",&ptimbalance_,"ptimbalance/F");
   tagprobe_tree_ -> Branch("dphi",&dphi_,"dphi/F");
   tagprobe_tree_ -> Branch("deta",&deta_,"deta/F");
   tagprobe_tree_ -> Branch("dr",&dr_,"dr/F");
   tagprobe_tree_ -> Branch("tag_qglikelihood",&tag_qglikelihood_,"tag_qglikelihood/F");
   tagprobe_tree_ -> Branch("tag_prob_b",&tag_prob_b_,"tag_prob_b/F");
   tagprobe_tree_ -> Branch("tag_prob_bb",&tag_prob_bb_,"tag_prob_bb/F");
   tagprobe_tree_ -> Branch("tag_prob_lepb",&tag_prob_lepb_,"tag_prob_lepb/F");
   tagprobe_tree_ -> Branch("tag_prob_c",&tag_prob_c_,"tag_prob_c_/F");
   tagprobe_tree_ -> Branch("tag_prob_g",&tag_prob_g_,"tag_prob_g/F");
   tagprobe_tree_ -> Branch("tag_prob_light",&tag_prob_light_,"tag_prob_light/F");
   tagprobe_tree_ -> Branch("tag_pt",&tag_pt_,"tag_pt/F");
   tagprobe_tree_ -> Branch("tag_eta",&tag_eta_,"tag_eta/F");
   tagprobe_tree_ -> Branch("probe_pt",&probe_pt_,"probe_pt/F");
   tagprobe_tree_ -> Branch("probe_eta",&probe_eta_,"probe_eta/F");
   tagprobe_tree_ -> Branch("weight",&tpweight_,"weight/F");
   
}

void OnlineBtagAnalyser::fillTagAndProbeTree()
{
   if ( ! do_tree_ ) return;
   std::string algo = config_->btagalgo_;
   if ( ! do_tree_ || algo != "deepflavour" ) return;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
      h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Fill TagAndProbe tree");
   
   int r1 = 1;
   int r2 = 2;
   int j1 = r1-1;
   int j2 = r2-1;
   
   if ( config_->isMC() && ( config_ -> useJetsFlavour() || config_ -> useJetsExtendedFlavour() ))
   {
      if ( config_ -> useJetsFlavour() )
      {
         probe_flavour_ = abs(selectedJets_[j1]->flavour());
         tag_flavour_   = abs(selectedJets_[j2]->flavour());
      }
      if ( config_ -> useJetsExtendedFlavour() )
      {
         std::string probe_flv = selectedJets_[j1]->extendedFlavour();
         std::string tag_flv   = selectedJets_[j2]->extendedFlavour();
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

   // variables for the tmva tree
   variables_ai(r1,r2);   
   
   tpweight_         = weight_;
   
   tagprobe_tree_ -> Fill();
   h1_["cutflow"] -> Fill(cutflow_,weight_);

}

void OnlineBtagAnalyser::variables_ai(const int & r1, const int & r2)
{
   int j1 = r1-1;
   int j2 = r2-1;
   ptimbalance_      = fabs(selectedJets_[j1]->pt() - selectedJets_[j2]->pt())/selectedJets_[j1]->pt();
   dphi_             = fabs(selectedJets_[j1]->deltaPhi(*selectedJets_[j2]));
   deta_             = fabs(selectedJets_[j1]->eta() - selectedJets_[j2]->eta());
   dr_               = fabs(selectedJets_[j1]->deltaR(*selectedJets_[j2]));
   tag_qglikelihood_ = selectedJets_[j2]->qgLikelihood();
   tag_prob_b_       = selectedJets_[j2]->btag("btag_dfb");
   tag_prob_bb_      = selectedJets_[j2]->btag("btag_dfbb");
   tag_prob_lepb_    = selectedJets_[j2]->btag("btag_dflepb");
   tag_prob_c_       = selectedJets_[j2]->btag("btag_dfc");
   tag_prob_g_       = selectedJets_[j2]->btag("btag_dfg");
   tag_prob_light_   = selectedJets_[j2]->btag("btag_dflight");
   tag_pt_           = selectedJets_[j2]->pt();
   tag_eta_          = selectedJets_[j2]->eta();
   probe_pt_         = selectedJets_[j1]->pt();
   probe_eta_        = selectedJets_[j1]->eta();
   
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
   variables_ai(1,2);
   
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
