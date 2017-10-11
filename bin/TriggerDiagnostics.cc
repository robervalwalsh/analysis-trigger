#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <boost/filesystem.hpp>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 
#include "TGraph.h" 
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TVectorD.h"

#include "Analysis/Core/interface/Analysis.h"

#include "HLTPathsAllHad.h"
#include "L1TriggersAllHad.h"
#include "L1TriggersSemiLep.h"

#include "lumis.h" 

using namespace std;
using namespace analysis;
using namespace analysis::core;


std::string inputList_;
std::string basename_;
bool isData_;

std::map<std::string, double> xsections_;
std::vector<std::string> jetTriggerObjects_;

// lumi / pile up from data
brilcalc bc_;
//std::map<int, std::map<int,float> > lumiByLS_;
std::map<int, std::map<int,float> > pileupByLS_;

// =============================================================================================   
int main(int argc, char * argv[])
{
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // Input files list
   inputList_ = "rootFileList.txt";
   isData_ = false;
   if ( argc == 2 )
   {
      inputList_ = std::string(argv[1]);
      basename_ =  std::string(boost::filesystem::basename(inputList_));
      std::cout << basename_ << std::endl;
      isData_ = (basename_.find("data") != std::string::npos);
   }
   
   Analysis analysis(inputList_,"MssmHbbTrigger/Events/EventInfo");
   
   // Physics Objects Collections
   // Trigger Results
   analysis.triggerResults("MssmHbbTrigger/Events/TriggerResults");
   
// Trigger objects
// L1 seeds
   jetTriggerObjects_.push_back("hltL1sSingleMu3");
   jetTriggerObjects_.push_back("hltL1sSingleJet20");
   
// HLT modules
   // L2 and L3 muons (HLT_Mu5)
   jetTriggerObjects_.push_back("hltL2fL1sMu3L1f0L2Filtered3Q");
   jetTriggerObjects_.push_back("hltL3fL1sMu3L1f0L2f3QL3Filtered5Q");
   
   // Calo and PF jets  (HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p014_DoublePFJet20)
   jetTriggerObjects_.push_back("hltBDiJet20L1FastJetCentral");
   jetTriggerObjects_.push_back("hltDoublePFBJet20");
   
   //
   
   std::string trgobj_path = "MssmHbbTrigger/Events/hltTriggerSummaryAOD/";
   for ( auto & obj : jetTriggerObjects_ )
      analysis.addTree<TriggerObject>(obj,trgobj_path+obj);
   
   // DATA only
   if ( isData_ )
   {
      bc_ = readLumisCsv("lumis.csv");
//      lumiByLS_ = bc_.lumiByLS;
      pileupByLS_ = bc_.pileupByLS;
      // Certified lumis
      analysis.processJsonFile("json.txt");
   }
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   int mainCounter = 0;
   int secCounter  = 0;
   int sec2Counter = 0;
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
      analysis.event(i);
      
      if ( i % 100000 == 0 && i>0 )
      {
         std::cout << i << " events analysed" << std::endl;
      }
      
      if ( analysis.triggerResult("HLT_BTagMu_DiJet20_Mu5_v") )
      {
         ++mainCounter;
         if ( analysis.triggerResult("HLT_Mu5_v") ) 
         {
            ++secCounter;
         }
         if ( analysis.triggerResult("HLT_L1SingleMu3_v") )
         {
            ++sec2Counter;
         }
      }
      
   }
   
   std::cout << "Main trigger accepts  = " << mainCounter << std::endl;
   std::cout << "Sec trigger accepts   = " << secCounter  << std::endl;
   std::cout << "Sec 2 trigger accepts = " << sec2Counter  << std::endl;
   std::cout << "Fraction   = " << double(secCounter)/mainCounter << std::endl;
   std::cout << "Fraction 2 = " << double(sec2Counter)/mainCounter << std::endl;
   
   
//    
}



// TRIGGER PATHS
//                                   'HLT_L1SingleJet20_v1',
//                                   'HLT_L1SingleMu3_v1',
//                                   'HLT_Mu5_v1',
//                                   'HLT_PFJet20_v1',
//                                   'HLT_PFJet40_v1',
//                                   'HLT_PFJet60_v1',
//                                   'HLT_DoublePFJet20_v1',
//                                   'HLT_DiPFJetAve40_v1',
//                                   'HLT_DiPFJetAve60_v1',
//                                   'HLT_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30_v1',
//                                   'HLT_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30MaxDeta1p6_v1',
//                                   'HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v1',
//                                   'HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_v1',
//                                   'HLT_DoubleJetsC100_SingleBTagCSV_p014_v1',
//                                   'HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_v1',
//                                   'HLT_BTagMu_DiJet20_Mu5_v1',
//                                   'HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p014_DoublePFJet20_v1',
//                                   'HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p014_DoublePFJet20MaxDeta1p6_v1',
//                                   'HLT_L1MuJet_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30_v1',
//                                   'HLT_L1MuJet_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30MaxDeta1p6_v1',



// HLT ConfDB table name = /users/rwalsh/dev/CMSSW_8_0_X/MssmHbb/V15
// GlobalTag = 80X_mcRun2_asymptotic_v20
// ========================================
// 
// L1RePack_step
// ----------------------------------------
// 
// HLTriggerFirstPath
// ----------------------------------------
// 
// HLT_L1SingleJet20_v1
//    L1T Seed: L1_SingleJet20,
//            Trigger Object:               'hltL1sSingleJet20',
// ----------------------------------------
// 
// HLT_L1SingleMu3_v1
//    L1T Seed: L1_SingleMu3,
//            Trigger Object:               'hltL1sSingleMu3',
//            Trigger Object:               'hltL1fL1sMu3L1Filtered0',
// ----------------------------------------
// 
// HLT_Mu5_v1
//    L1T Seed: L1_SingleMu3,
//            Trigger Object:               'hltL1sSingleMu3',
//            Trigger Object:               'hltL1fL1sMu3L1Filtered0',
//            Trigger Object:               'hltL2fL1sMu3L1f0L2Filtered3Q',
//            Trigger Object:               'hltL3fL1sMu3L1f0L2f3QL3Filtered5Q',
// ----------------------------------------
// 
// HLT_PFJet20_v1
//    L1T Seed: L1_SingleJet20,
//            Trigger Object:               'hltL1sSingleJet20',
//            Trigger Object:               'hltSingleCaloJet20',
//            Trigger Object:               'hltSinglePFJet20',
// ----------------------------------------
// 
// HLT_PFJet40_v1
//    L1T Seed: L1_ZeroBias,
//            Trigger Object:               'hltL1sZeroBias',
//            Trigger Object:               'hltSingleCaloJet10',
//            Trigger Object:               'hltSinglePFJet40',
// ----------------------------------------
// 
// HLT_PFJet60_v1
//    L1T Seed: L1_SingleJet35,
//            Trigger Object:               'hltL1sSingleJet35',
//            Trigger Object:               'hltSingleCaloJet40',
//            Trigger Object:               'hltSinglePFJet60',
// ----------------------------------------
// 
// HLT_DoublePFJet20_v1
//    L1T Seed: L1_SingleJet20,
//            Trigger Object:               'hltL1sSingleJet20',
//            Trigger Object:               'hltDoubleCaloJet20',
//            Trigger Object:               'hltDoublePFJet20',
// ----------------------------------------
// 
// HLT_DiPFJetAve40_v1
//    L1T Seed: L1_ZeroBias,
//            Trigger Object:               'hltL1sZeroBias',
//            Trigger Object:               'hltDiCaloJetAve30',
//            Trigger Object:               'hltDiPFJetAve40',
// ----------------------------------------
// 
// HLT_DiPFJetAve60_v1
//    L1T Seed: L1_ZeroBias,
//            Trigger Object:               'hltL1sZeroBias',
//            Trigger Object:               'hltDiCaloJetAve45',
//            Trigger Object:               'hltDiPFJetAve60',
// ----------------------------------------
// 
// HLT_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30_v1
//    L1T Seed: L1_SingleJet20,
//            Trigger Object:               'hltL1sSingleJet20',
//            Trigger Object:               'hltDoubleCaloJet30',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJet30',
// ----------------------------------------
// 
// HLT_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30MaxDeta1p6_v1
//    L1T Seed: L1_SingleJet20,
//            Trigger Object:               'hltL1sSingleJet20',
//            Trigger Object:               'hltDoubleCaloJet30',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJet30',
//            Trigger Object:               'hltDoublePFJetsC30MaxDeta1p6',
// ----------------------------------------
// 
// HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v1
//    L1T Seed: L1_DoubleJetC100,
//            Trigger Object:               'hltL1sDoubleJetC100',
//            Trigger Object:               'hltDoubleJetsC100',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJetsC100',
//            Trigger Object:               'hltDoublePFJetsC100MaxDeta1p6',
// ----------------------------------------
// 
// HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6_v1
//    L1T Seed: L1_DoubleJetC112,
//            Trigger Object:               'hltL1sDoubleJetC112',
//            Trigger Object:               'hltDoubleJetsC112',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJetsC112',
//            Trigger Object:               'hltDoublePFJetsC112MaxDeta1p6',
// ----------------------------------------
// 
// HLT_DoubleJetsC100_SingleBTagCSV_p014_v1
//    L1T Seed: L1_DoubleJetC100,
//            Trigger Object:               'hltL1sDoubleJetC100',
//            Trigger Object:               'hltDoubleJetsC100',
//            Trigger Object:               'hltSingleBTagCSV0p84',
// ----------------------------------------
// 
// HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350_v1
//    L1T Seed: L1_DoubleJetC100,
//            Trigger Object:               'hltL1sDoubleJetC100',
//            Trigger Object:               'hltDoubleJetsC100',
//            Trigger Object:               'hltSingleBTagCSV0p84',
//            Trigger Object:               'hltJetC350',
// ----------------------------------------
// 
// HLT_BTagMu_DiJet20_Mu5_v1
//    L1T Seed: L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4,
//            Trigger Object:               'hltL1sMu3JetC16dEtaMax0p4dPhiMax0p4',
//            Trigger Object:               'hltBDiJet20L1FastJetCentral',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetL25FilterByDR',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetMu5L3FilterByDR',
// ----------------------------------------
// 
// HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p014_DoublePFJet20_v1
//    L1T Seed: L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4,
//            Trigger Object:               'hltL1sMu3JetC16dEtaMax0p4dPhiMax0p4',
//            Trigger Object:               'hltBDiJet20L1FastJetCentral',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetL25FilterByDR',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetMu5L3FilterByDR',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFBJet20',
// ----------------------------------------
// 
// HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p014_DoublePFJet20MaxDeta1p6_v1
//    L1T Seed: L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4,
//            Trigger Object:               'hltL1sMu3JetC16dEtaMax0p4dPhiMax0p4',
//            Trigger Object:               'hltBDiJet20L1FastJetCentral',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetL25FilterByDR',
//            Trigger Object:               'hltBSoftMuonDiJet20L1FastJetMu5L3FilterByDR',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFBJet20',
//            Trigger Object:               'hltDoublePFJetsC20MaxDeta1p6',
// ----------------------------------------
// 
// HLT_L1MuJet_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30_v1
//    L1T Seed: L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4,
//            Trigger Object:               'hltL1sMu3JetC16dEtaMax0p4dPhiMax0p4',
//            Trigger Object:               'hltDoubleCaloJet30',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJet30',
// ----------------------------------------
// 
// HLT_L1MuJet_DoubleJet30_DoubleBTagCSV_p014_DoublePFJet30MaxDeta1p6_v1
//    L1T Seed: L1_Mu3_JetC16_dEta_Max0p4_dPhi_Max0p4,
//            Trigger Object:               'hltL1sMu3JetC16dEtaMax0p4dPhiMax0p4',
//            Trigger Object:               'hltDoubleCaloJet30',
//            Trigger Object:               'hltBTagCaloCSVp014DoubleWithMatching',
//            Trigger Object:               'hltDoublePFJet30',
//            Trigger Object:               'hltDoublePFJetsC30MaxDeta1p6',
// ----------------------------------------
// 
// HLTriggerFinalPath
// ----------------------------------------
// 
