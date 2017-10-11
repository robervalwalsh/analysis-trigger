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
#include "TLorentzVector.h"

#include "Analysis/Core/interface/Analysis.h"

#include "lumis.h" 

using namespace std;
using namespace analysis;
using namespace analysis::core;


std::map<std::string, TH1F*> h1_;
std::map<std::string, TH1F*> h1x_;

TGraphAsymmErrors * g_eff_;
TGraphErrors * g_rates_;
std::string inputList_;
std::string basename_;
std::vector<std::string> triggers_;
std::map<std::string, bool> trigger_accept_;

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

double etacut_ = 2.3;
double detacut_ = 1.6;

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   xsections_ = CrossSections();
   
   isData_ = false;
   // Input files list
   inputList_ = "rootFileList.txt";
   if ( argc == 2 )
   {
      inputList_ = std::string(argv[1]);
   }
   
   basename_ =  std::string(boost::filesystem::basename(inputList_));
   std::cout << basename_ << std::endl;
   TNamed sampleName("SampleName",basename_.c_str());
   std::string basename_orginal = basename_;
   isData_ = (basename_.find("data") != std::string::npos);
   if ( (basename_.find("_x") != std::string::npos) )
   {
      basename_.erase(basename_.length()-5,5);
   }
   
   std::string  refTrigger = "";
   if ( isData_ )
      refTrigger = "HLT_CaloJets_Muons_CaloBTagCSV_PFJets_v";
   else
      refTrigger = "HLT_ZeroBias_v";
   
   if ( isData_ ) std::cout << "Running on Data"  << std::endl;
   else           std::cout << "Running on Simulation"  << std::endl;
   
   int ncoll = 2200;
   int ncollFix = -1;

   xsections_ = CrossSections();
   
   triggers_.push_back("ZeroBias");
   
   
   // semileptonic
   trigger_accept_["L1_Mu12_DoubleJet40"] = true;
   trigger_accept_["L1_Mu12_DoubleJet40dEta1p6"] = true;
   trigger_accept_["L1_Mu3_Jet16"] = true;
   trigger_accept_["L1_Mu3_Jet60"] = true;
   trigger_accept_["L1_Mu3_Jet120"] = true;
   trigger_accept_["L1_DoubleJet40"] = true;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] = false;
   // backup
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"] = false;
   
   // monitor
   trigger_accept_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] = false;
   
   
   // alhadronic
   trigger_accept_["L1_DoubjeJet100"] = true;
   trigger_accept_["L1_DoubjeJet112"] = true;
   trigger_accept_["L1_DoubjeJet100dEta1p6"] = true;
   trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"] = false;
   trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"] = false;
   // backup
   trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110"] = false;
   trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120"] = false;
   trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130"] = false;
   
     // monitor
   trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100"] = false;
   trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140"] = false;
   trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200"] = false;
   trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350"] = false;
   
   for ( auto & ta : trigger_accept_ )
   {
      if ( ta.first.find("HLT") != std::string::npos )  triggers_.push_back(ta.first);
   }
   
   for ( size_t i = 0 ; i < triggers_.size(); ++i )
   {
      h1_[triggers_[i]]                        = new TH1F(Form("h_n%s",triggers_[i].c_str())    ,"",70,0.,70.);
      h1_[Form("%s_psw",triggers_[i].c_str())] = new TH1F(Form("h_n%s_psw",triggers_[i].c_str()),"",70,0.,70.);
   }   
   
   h1_["Lumis"] = new TH1F("Lumis","lumi sections per pile-up",70,0.,70.);
   
   
   Analysis analysis(inputList_,"MssmHbbTrigger/Events/EventInfo");
   
   // Physics Objects Collections
   // Trigger
   analysis.triggerResults("MssmHbbTrigger/Events/TriggerResults");
   
   // Get cross section
   double crossSection  = -1.;
   if ( ! isData_ )
   {
      analysis.crossSections("MssmHbbTrigger/Metadata/CrossSections");
      crossSection = analysis.crossSection();
      crossSection = xsections_[basename_];
   }
   
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
      
      int run = analysis.run();
      int ls  = analysis.lumiSection();
      
      float nPileup;
      if ( isData_ )
      {
         if ( ! analysis.selectJson() ) continue;
         nPileup = pileupByLS_[run][ls];
         prescale = lumiByLS_[run][ls]/lumiPathByLS_[run][ls];
         ncoll = coll_[run][ls];
         
         // {"283171": [[56,69],[126,144],[175,193],[209,233]]}
         if ( run == 283171 )
         {
            prescale = 3606;                                                                           
            if ( ls >=  56 && ls <=  69 ) nPileup = 55;                                                
            if ( ls >= 126 && ls <= 144 ) nPileup = 50;                                                
            if ( ls >= 175 && ls <= 193 ) nPileup = 47;                                                
            if ( ls >= 209 && ls <= 233 ) nPileup = 45;                                                
         }                                                                                             
                                                                                                       
         if ( ncollFix > 0 && ncoll != ncollFix ) continue;  // I think that's for tests only, check   
                                                                                                       
         // Lumi sections counter (for data)                                                           
         if ( ! lsXpu_[run][ls] )                                                                      
         {                                                                                             
            lsXpu_[run][ls] = true;                                                                    
            h1_["Lumis"] -> Fill(nPileup);                                                             
         }                                                                                             
      }                                                                                                
      else                                                                                             
      {                                                                                                
         nPileup = (float)analysis.nTruePileup(); // from eventinfo MC                                 
         ncoll = 1;                                                                                    
      }                                                                                                
                                                                                                       
      // ZEROBIAS                                                                                      
      h1_["ZeroBias"] -> Fill(nPileup, 1./ncoll);                                                      
      h1_["ZeroBias_psw"] -> Fill(nPileup, prescale/ncoll);                                            
                                                                                                       
                                                                                                       
// ======== Begin of L1 ======== //      

// ======== End of L1 ======== //      
      
// ======== Begin of HLT ======== //  

      // All hadronic                                                                                       
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"]       =    false;                      
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"]       =    analysis.triggerResult("HLT_DoubleJets100_DoubleBTagCSV_0p92_DoublePFJets100MaxDeta1p6_v");
      // monitor                                                                                            
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100"]       =    analysis.triggerResult("HLT_DoubleJets100_SingleBTagCSV_0p92_DoublePFJets100_v");
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140"]       =    false;                      
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200"]       =    analysis.triggerResult("HLT_DoubleJets100_SingleBTagCSV_0p92_DoublePFJets200_v");
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350"]       =    analysis.triggerResult("HLT_DoubleJets100_SingleBTagCSV_0p92_DoublePFJets350_v");
                                                                                                            
      // backup                                                                                             
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110"]       =    analysis.triggerResult("HLT_DoubleJets100_DoubleBTagCSV_0p92_DoublePFJets116MaxDeta1p6_v");
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120"]       =    false;
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130"]       =    analysis.triggerResult("HLT_DoubleJets100_DoubleBTagCSV_0p92_DoublePFJets128MaxDeta1p6_v");
                                                                                                            
      // Semileptonic                                                                                       
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"]  =    false;
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"]  =    analysis.triggerResult("HLT_DoubleJets30_Mu12_DoubleBTagCSV_0p92_DoublePFJets40MaxDeta1p6_v");
      // monitor                                                                                            
      trigger_accept_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"]         =    analysis.triggerResult("HLT_DoubleJets30_SingleBTagCSV_0p92_DoublePFJets40_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"]  =    analysis.triggerResult("HLT_DoubleJets30_Mu12_SingleBTagCSV_0p92_DoublePFJets40_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] =    analysis.triggerResult("HLT_DoubleJets30_Mu12_SingleBTagCSV_0p92_DoublePFJets100_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] =    analysis.triggerResult("HLT_DoubleJets30_Mu12_SingleBTagCSV_0p92_DoublePFJets200_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] =    analysis.triggerResult("HLT_DoubleJets30_Mu12_SingleBTagCSV_0p92_DoublePFJets350_v");
      // backup                                                                                             
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"]  =    analysis.triggerResult("HLT_DoubleJets30_Mu12_DoubleBTagCSV_0p92_DoublePFJets54MaxDeta1p6_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"]  =    analysis.triggerResult("HLT_DoubleJets30_Mu12_DoubleBTagCSV_0p92_DoublePFJets62MaxDeta1p6_v");
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"]  =    false;
      
      
// --------------------      
      // HLT Path - All had
      if ( trigger_accept_["L1_DoubjeJet100dEta1p6"] )
      {
         if ( trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"] )
         {
            h1_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
            h1_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"] )
         {
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         // backup
         if ( trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110"] )
         {
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110"] -> Fill(nPileup, 1./ncoll);   
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120"] )
         {
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120"] -> Fill(nPileup, 1./ncoll);   
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130"] )
         {
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130"] -> Fill(nPileup, 1./ncoll);   
            h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         
      }
      // All had monitor
      if ( trigger_accept_["L1_DoubjeJet100"] && trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100"] )
      {
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_DoubjeJet112"] && trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140"] )
      {
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_DoubjeJet112"] && trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200"] )
      {
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_DoubjeJet112"] && trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350"] )
      {
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      
            
      // HLT Path - Semileptonic
      if ( trigger_accept_["L1_Mu12_DoubleJet40dEta1p6"] )
      {
         if ( trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] )
         {
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] -> Fill(nPileup, 1./ncoll);
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] )
         {
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] -> Fill(nPileup, 1./ncoll);
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         // backup
         if ( trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"] )
         {
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"] -> Fill(nPileup, 1./ncoll);
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"] )
         {
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"] -> Fill(nPileup, 1./ncoll);
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
         if ( trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"] )
         {
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"] -> Fill(nPileup, 1./ncoll);
            h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70_psw"] -> Fill(nPileup, prescale/ncoll);   
         }
      }
      // semilep monitor
      if ( trigger_accept_["L1_DoubleJet40"] && trigger_accept_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"] )
      {
         h1_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_Mu3_Jet16"] && trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"] )
      {
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_Mu3_Jet60"] && trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] )
      {
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_Mu3_Jet120"] && trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] )
      {
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( trigger_accept_["L1_Mu3_Jet120"] && trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] )
      {
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      
// ======== End of HLT ======== //  
      
   }
   
//    // number of zerobias events
   TVectorD nZB(1);
   if ( ! isData_ )
      nZB[0] = h1_["ZeroBias"] -> GetEntries();
   else
      nZB[0] = -1;
//    
   TFile * f_out = new TFile(Form("mssmhbb_triggers_%s.root",basename_orginal.c_str()),"RECREATE");
   for ( auto & h : h1_ )
   {
      if ( ! isData_ && h.first == "Lumis" ) continue;
      std::cout << h.first << std::endl;
      h.second->Write();
   }
   
   TVectorD xsection(1);
   xsection[0] = crossSection;
   std::cout << "Cross section = " << crossSection << std::endl;
   xsection.Write("xsection");
   
   sampleName.Write();
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
