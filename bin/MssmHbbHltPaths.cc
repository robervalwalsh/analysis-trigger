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
using namespace analysis::tools;


std::map<std::string, TH1F*> h1_;
std::map<std::string, TH1F*> h1x_;

TGraphAsymmErrors * g_eff_;
TGraphErrors * g_rates_;
std::string inputList_;
std::string basename_;
std::vector<std::string> triggers_;
bool isData_;

std::map<std::string, double> xsections_;
std::vector<std::string> jetTriggerObjects_;

// lumi / pile up from data
brilcalc bc_;
brilcalc bcp_;
collisions coll_;
std::map<int, std::map<int,float> > lumiByLS_;
std::map<int, std::map<int,float> > pileupByLS_;
std::map<int, std::map<int,float> > lumiPathByLS_;

std::map<int, std::map<int,bool> > lsXpu_;

std::map<std::string, double> CrossSections();

// =============================================================================================   
int main(int argc, char * argv[])
{
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   std::string  refTrigger = "";
//   refTrigger = "HLT_BTagMu_DiJet20_Mu5_v";
   refTrigger = "HLT_BTagMu_DiJet20_Mu5_DoubleBTagCSV_p000_DoublePFJet20_v";
   
   int ncoll = 2200;
   int ncollFix = -1;

   xsections_ = CrossSections();
   
   triggers_.push_back("BTagMu20");
   triggers_.push_back("BTagMu20BTagCSV000");
   triggers_.push_back("BTagMuDijet40Mu12BTagCSV000Mu5NewL1");
   triggers_.push_back("BTagMuDijet40Mu12BTagCSV084Mu5NewL1");
   triggers_.push_back("BTagMuDijet40Mu12BTagCSV088Mu5NewL1");
   triggers_.push_back("BTagMuDijet40Mu12BTagCSV092Mu5NewL1");
   triggers_.push_back("BTagMuDijet40Mu12BTagCSV094Mu5NewL1");
   triggers_.push_back("BTagMuDijet50&40Mu12BTagCSV084Mu5NewL1");
   
   for ( size_t i = 0 ; i < triggers_.size(); ++i )
   {
      h1_[triggers_[i]]                        = new TH1F(Form("h_n%s",triggers_[i].c_str())    ,"",70,0.,70.);
      h1_[Form("%s_psw",triggers_[i].c_str())] = new TH1F(Form("h_n%s_psw",triggers_[i].c_str()),"",70,0.,70.);
      // distributions
      h1x_[Form("jetsN_%s",triggers_[i].c_str())]   = new TH1F(Form("h_jetsN_%s",triggers_[i].c_str()),"",10,0,10);
      h1x_[Form("jetsPT_%s",triggers_[i].c_str())]  = new TH1F(Form("h_jetsPT_%s",triggers_[i].c_str()),"",2000,0,1000);
      h1x_[Form("jetsETA_%s",triggers_[i].c_str())] = new TH1F(Form("h_jetsETA_%s",triggers_[i].c_str()),"",500,-5,5);
      h1x_[Form("jetsPHI_%s",triggers_[i].c_str())] = new TH1F(Form("h_jetsPHI_%s",triggers_[i].c_str()),"",630,-3.15,3.15);
      h1x_[Form("jet1PT_%s",triggers_[i].c_str())]  = new TH1F(Form("h_jet1PT_%s",triggers_[i].c_str()),"",2000,0,1000);
      h1x_[Form("jet1ETA_%s",triggers_[i].c_str())] = new TH1F(Form("h_jet1ETA_%s",triggers_[i].c_str()),"",500,-5,5);
      h1x_[Form("jet1PHI_%s",triggers_[i].c_str())] = new TH1F(Form("h_jet1PHI_%s",triggers_[i].c_str()),"",630,-3.15,3.15);

   }
   
   h1_["Lumis"] = new TH1F("Lumis","lumi sections per pile-up",70,0.,70.);
   
   // Input files list
   inputList_ = "rootFileList.txt";
   isData_ = true;
   if ( argc == 2 )
   {
      inputList_ = std::string(argv[1]);
      basename_ =  std::string(boost::filesystem::basename(inputList_));
      std::cout << basename_ << std::endl;
      isData_ = (basename_.find("data") != std::string::npos);
   }
   
   Analysis analysis(inputList_,"MssmHbbTrigger/Events/EventInfo");
   
   // Physics Objects Collections
   analysis.triggerResults("MssmHbbTrigger/Events/TriggerResults");
   if ( refTrigger != "HLT_BTagMu_DiJet20_Mu5_v" )
      analysis.addTree<JetTag> ("JetsTags","MssmHbbTrigger/Events/hltCombinedSecondaryVertexBJetTagsCalo");
   
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
   jetTriggerObjects_.push_back("hltDoublePFJet20");
   
   //
   
   std::string trgobj_path = "MssmHbbTrigger/Events/hltTriggerSummaryAOD/";
   for ( auto & obj : jetTriggerObjects_ )
      analysis.addTree<TriggerObject>(obj,trgobj_path+obj);
   
   // Get cross section
//    double crossSection  = -1.;
//    if ( ! isData_ )
//    {
//       analysis.crossSections("MssmHbbTrigger/Metadata/CrossSections");
//       crossSection = analysis.crossSection();
//       crossSection = xsections_[basename_];
//    }
   
   // DATA only
   if ( isData_ )
   {
      bc_ = readLumisCsv("lumis.csv");
      bcp_ = readLumisPathCsv("lumispath.csv");
      coll_ = readCollisionsCsv("collisions.csv");
      pileupByLS_ = bc_.pileupByLS;
      lumiByLS_ = bc_.lumiByLS;
      lumiPathByLS_ = bcp_.lumiByLS;
      // Lumi sections versus pile-up
      for ( auto & i1 : pileupByLS_ )
      {
         // int run = i1.first;
         for ( auto & i2 : i1.second )
         {
            // int lumisection = i2.first;
            lsXpu_[i1.first][i2.first] = false;
         }
      }
      // Certified lumis
      analysis.processJsonFile("json.txt");
   }
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
      analysis.event(i);
      
      if ( i % 100000 == 0 && i>0 )
      {
         std::cout << i << " events analysed" << std::endl;
      }
      
      // hltPath0 - reference trigger
      if ( ! analysis.triggerResult(refTrigger) ) continue;
      
      float prescale = 1.;
      
      float nPileup;
      if ( isData_ )
      {
         if ( ! analysis.selectJson() ) continue;
         nPileup = pileupByLS_[analysis.run()][analysis.lumiSection()];
         prescale = lumiByLS_[analysis.run()][analysis.lumiSection()]/lumiPathByLS_[analysis.run()][analysis.lumiSection()];
         ncoll = coll_[analysis.run()][analysis.lumiSection()];
//         std::cout << analysis.run() << "   " << analysis.lumiSection() << "   " << ncoll << std::endl;
         if ( ncollFix > 0 && ncoll != ncollFix ) continue;
         if ( ! lsXpu_[analysis.run()][analysis.lumiSection()] )
         {
//            std::cout << analysis.run() << " , " << analysis.lumiSection() << "  ps = " << prescale << " pile up = " << nPileup <<  std::endl;
            lsXpu_[analysis.run()][analysis.lumiSection()] = true;
            h1_["Lumis"] -> Fill(nPileup);
         }
      }
      else
      {
         nPileup = (float)analysis.nTruePileup();
         ncoll = 1;
      }
      
      std::shared_ptr< Collection<JetTag> > jetstags;
      
      if ( refTrigger != "HLT_BTagMu_DiJet20_Mu5_v" ) 
      {
         // btagging - up to 6 jets with pt > 30 GeV
         jetstags = analysis.collection<JetTag>("JetsTags");
         if ( analysis.triggerResult(refTrigger) && ncoll > 0 )
         {
            if (  L1Mu10DiJet32(analysis) && analysis.triggerResult("HLT_Mu5_v") )
            {
               // select muons with pt > 10 eta < 2.3
               auto hlt_mu5 = analysis.collection<TriggerObject>("hltL3fL1sMu3L1f0L2f3QL3Filtered5Q");
               std::vector<TriggerObject> hltmu12;
               for ( int m = 0; m < hlt_mu5->size() ; ++m )
               {
                  TriggerObject mu = hlt_mu5->at(m);
                  if ( mu.pt() < 12. || fabs(mu.eta()) > 2.3 ) continue;
                  hltmu12.push_back(mu);
               }
               auto hlt_2pfjet20 = analysis.collection<TriggerObject>("hltDoublePFBJet20");
               std::vector<TriggerObject> hltpfjet40;
               for ( int j = 0; j < hlt_2pfjet20->size() ; ++j )
               {
                  TriggerObject jet = hlt_2pfjet20->at(j);
                  if ( jet.pt() < 40. || fabs(jet.eta()) > 2.3 ) continue;
                  hltpfjet40.push_back(jet);
               }
               std::vector<TriggerObject> hltpfjet50;
               for ( int j = 0; j < hlt_2pfjet20->size() ; ++j )
               {
                  TriggerObject jet = hlt_2pfjet20->at(j);
                  if ( jet.pt() < 50. || fabs(jet.eta()) > 2.3 ) continue;
                  hltpfjet50.push_back(jet);
               }
               
               
               std::vector<JetTag> btag084jets;
               std::vector<JetTag> btag088jets;
               std::vector<JetTag> btag092jets;
               std::vector<JetTag> btag094jets;
               
//               std::cout <<  "oioi   " << jetstags << "   "  << jetstags->size() << std::endl;
               for ( int j = 0; j < jetstags->size() ; ++j )
               {
                  if ( j > 5 ) break;
                  JetTag jet = jetstags->at(j);
                  if (jet.btag() >= 0.84 )  btag084jets.push_back(jet);
                  if (jet.btag() >= 0.88 )  btag088jets.push_back(jet);
                  if (jet.btag() >= 0.92 )  btag092jets.push_back(jet);
                  if (jet.btag() >= 0.94 )  btag094jets.push_back(jet);
               }
               if ( hltmu12.size() >= 1 && hltpfjet40.size() >= 2 )
               {
                  h1_["BTagMuDijet40Mu12BTagCSV000Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                  h1_["BTagMuDijet40Mu12BTagCSV000Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                  if ( btag084jets.size() >= 2 )
                  {
                     h1_["BTagMuDijet40Mu12BTagCSV084Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                     h1_["BTagMuDijet40Mu12BTagCSV084Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                     if ( hltpfjet50.size() >= 1 )
                     {
                       h1_["BTagMuDijet50&40Mu12BTagCSV084Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                       h1_["BTagMuDijet50&40Mu12BTagCSV084Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                     }
                  }
                  if ( btag088jets.size() >= 2 )
                  {
                     h1_["BTagMuDijet40Mu12BTagCSV088Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                     h1_["BTagMuDijet40Mu12BTagCSV088Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                  }
                  if ( btag092jets.size() >= 2 )
                  {
                     h1_["BTagMuDijet40Mu12BTagCSV092Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                     h1_["BTagMuDijet40Mu12BTagCSV092Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                  }
                  if ( btag094jets.size() >= 2 )
                  {
                     h1_["BTagMuDijet40Mu12BTagCSV094Mu5NewL1"] -> Fill(nPileup, 1./ncoll);
                     h1_["BTagMuDijet40Mu12BTagCSV094Mu5NewL1_psw"] -> Fill(nPileup, prescale/ncoll);
                  }
               }
            }
         }
         
      }
      
      // hltPath0 - reference trigger
      if ( analysis.triggerResult("HLT_BTagMu_DiJet20_Mu5_v") && ncoll > 0 ) 
      {
         h1_["BTagMu20"] -> Fill(nPileup, 1./ncoll);
         h1_["BTagMu20_psw"] -> Fill(nPileup, prescale/ncoll);
      }
      
      if ( analysis.triggerResult(refTrigger) && ncoll > 0 ) 
      {
         h1_["BTagMu20BTagCSV000"] -> Fill(nPileup, 1./ncoll);
         h1_["BTagMu20BTagCSV000_psw"] -> Fill(nPileup, prescale/ncoll);
      }
   }
   
//    // number of zerobias events
//    TVectorD nZB(1);
//    nZB[0] = h1_["ZeroBias"] -> GetEntries();
//    
   TFile * f_out = new TFile(Form("mssmhbb_triggers_%s.root",basename_.c_str()),"RECREATE");
   for ( auto & h : h1_ )
   {
      if ( ! isData_ && h.first == "Lumis" ) continue;
      std::cout << h.first << std::endl;
      h.second->Write();
   }
//    
//    TVectorD xsection(1);
//    xsection[0] = crossSection;
//    xsection.Write("xsection");
//    
//    TNamed sampleName("SampleName",basename_.c_str());
//    sampleName.Write();
//    
   f_out -> Close();
//    
//    TFile * f_dist = new TFile(Form("mssmhbb_triggers_distributions_%s.root",basename_.c_str()),"RECREATE");
//    
//    for ( auto & h : h1x_ )
//    {
//       std::cout << h.first << std::endl;
//       h.second->Write();
//    }
// 
//    xsection.Write("xsection");
//    nZB.Write("nZeroBias");
//    
//    f_dist -> Close();   
   
//    
}


std::map<std::string, double> CrossSections()
{
   std::map<std::string, double> xsections;
   xsections["QCD_Pt_15to30"]   = 1.83741e+09;
   xsections["QCD_Pt_30to50"]   = 1.40932e+08;
   xsections["QCD_Pt_50to80"]   = 1.92043e+07;
   xsections["QCD_Pt_80to120"]  = 2.76253e+06;
   xsections["QCD_Pt_120to170"] = 471100;
   xsections["QCD_Pt_170to300"] = 117276;
   xsections["QCD_Pt_300to470"] = 7823.28;
   xsections["QCD_Pt_470to600"] = 648.174;
   return xsections;
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
