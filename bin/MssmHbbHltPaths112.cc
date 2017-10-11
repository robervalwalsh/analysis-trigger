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
   trigger_accept_["L1_Mu12_DoubleJet40"] = false;
   trigger_accept_["L1_Mu12_DoubleJet40dEta1p6"] = false;
   trigger_accept_["L1_Mu3_Jet16"] = false;
   trigger_accept_["L1_Mu3_Jet60"] = false;
   trigger_accept_["L1_Mu3_Jet120"] = false;
   trigger_accept_["L1_DoubleJet40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] = false;
   // backup
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"] = false;
   
   // monitor
   trigger_accept_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"] = false;
   trigger_accept_["HLT_1CaloJets30_Muon12_1PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] = false;
   trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] = false;
   
   
   // alhadronic
   trigger_accept_["L1_DoubjeJet100"] = false;
   trigger_accept_["L1_DoubjeJet112"] = false;
   trigger_accept_["L1_DoubjeJet100dEta1p6"] = false;
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
   // Vertices
   analysis.addTree<Vertex> ("FastPV","MssmHbbTrigger/Events/hltFastPrimaryVertex");
   analysis.addTree<Vertex> ("FastPVPixel","MssmHbbTrigger/Events/hltFastPVPixelVertices");
   
   // Jets
   analysis.addTree<Jet> ("L1Jets","MssmHbbTrigger/Events/l1tJets");
   analysis.addTree<Jet> ("CaloJets","MssmHbbTrigger/Events/hltAK4CaloJetsCorrectedIDPassed");
   analysis.addTree<Jet> ("PFJets","MssmHbbTrigger/Events/hltAK4PFJetsLooseIDCorrected");
   
   // Muons
   analysis.addTree<Muon> ("L1Muons","MssmHbbTrigger/Events/l1tMuons");
   analysis.addTree<Muon> ("L2Muons","MssmHbbTrigger/Events/hltL2MuonCandidates");
   analysis.addTree<Muon> ("L3Muons","MssmHbbTrigger/Events/hltL3MuonCandidates");
   
   // BTag
   analysis.addTree<JetTag> ("JetsTags","MssmHbbTrigger/Events/hltCombinedSecondaryVertexBJetTagsCalo");

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
      
      
      // primary vertex
      std::shared_ptr< Collection<Vertex> > fastPVs = analysis.collection<Vertex>("FastPV");
      std::shared_ptr< Collection<Vertex> > fastPVsPixel = analysis.collection<Vertex>("FastPVPixel");
      if ( fastPVs -> size() < 1 || fastPVsPixel -> size() < 1 ) continue;
      Vertex fastpv = fastPVs->at(0);
      if ( !( ! fastpv.fake() && fastpv.ndof() > 0 && fabs(fastpv.z()) <= 25 && fastpv.rho() <= 2 ) ) continue;
      Vertex fastpvpix = fastPVsPixel->at(0);
      if ( !( ! fastpvpix.fake() && fastpvpix.ndof() > 0 && fabs(fastpvpix.z()) <= 25 && fastpvpix.rho() <= 2 ) ) continue;
      
// ======== Begin of L1 ======== //      

// --------------------      
// L1 Jets collection
      std::shared_ptr< Collection<Jet> > l1jets = analysis.collection<Jet>("L1Jets");
// L1 Muons collection
      std::shared_ptr< Collection<Muon> > l1muons = analysis.collection<Muon>("L1Muons");
// --------------------      
// L1 jets objects
      std::vector<Jet> l1jet16;
      std::vector<Jet> l1jet40;
      std::vector<Jet> l1jet60;
      std::vector<Jet> l1jet100;
      std::vector<Jet> l1jet112;
      std::vector<Jet> l1jet120;
// L1 delta_eta jets objects
      std::vector<Jet> l1jet40deta;
      std::vector<Jet> l1jet100deta;
// L1 Muons objects
      std::vector<Muon> l1mu3;
      std::vector<Muon> l1mu12;
// L1 jet-mu objects
      std::vector<Jet> l1muonjets;
      std::vector<Muon> l1jetmuons;
      std::vector<Jet> l1muonjetsbtv16;
      std::vector<Muon> l1jetmuonsbtv16;
      std::vector<Jet> l1muonjetsbtv60;
      std::vector<Muon> l1jetmuonsbtv60;
      std::vector<Jet> l1muonjetsbtv120;
      std::vector<Muon> l1jetmuonsbtv120;
// --------------------      
      // L1 jets selection
      for ( int j = 0; j < l1jets->size() ; ++j )
      {
         Jet l1j = l1jets->at(j);
         if ( l1j.pt() >= 16  && fabs(l1j.eta()) <= 3.0 ) l1jet16.push_back(l1j);
         if ( l1j.pt() >= 40  && fabs(l1j.eta()) <= etacut_ ) l1jet40.push_back(l1j);
         if ( l1j.pt() >= 60  && fabs(l1j.eta()) <= etacut_ ) l1jet60.push_back(l1j);
         if ( l1j.pt() >= 100 && fabs(l1j.eta()) <= etacut_ ) l1jet100.push_back(l1j);
         if ( l1j.pt() >= 112 && fabs(l1j.eta()) <= etacut_ ) l1jet112.push_back(l1j);
         if ( l1j.pt() >= 120 && fabs(l1j.eta()) <= etacut_ ) l1jet120.push_back(l1j);
      }
      // L1 delta_eta jets selection
      if ( l1jet40.size() >= 2 )
      {
         for ( int j1 = 0 ; j1 < int(l1jet40.size())-1 ; ++j1 )
         {
            Jet jet1 = l1jet40.at(j1);
            for ( int j2 = j1+1 ; j2 < int(l1jet40.size()) ; ++j2 )
            {
               Jet jet2 = l1jet40.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= detacut_ )
               {
                  l1jet40deta.push_back(jet1);
                  l1jet40deta.push_back(jet2);
                  if ( jet1.pt() >= 100 && jet2.pt() >= 100 ) // avoid repetition
                  {
                     l1jet100deta.push_back(jet1);
                     l1jet100deta.push_back(jet2);
                  }
               }
            }
         }
      }
      // L1 muon selection
      for ( int m = 0; m < l1muons->size() ; ++m )
      {
         Muon l1m = l1muons->at(m);
         if ( l1m.pt() >= 12  && fabs(l1m.eta()) <= etacut_ ) l1mu12.push_back(l1m);
         if ( l1m.pt() >= 3  && fabs(l1m.eta()) <= 2.5 ) l1mu3.push_back(l1m);
      }
      // L1 muon-jet selection
      if ( l1jet40.size() >= 1 && l1mu12.size() >= 1 )
      {
         for ( int j = 0 ; j < int(l1jet40.size()) ; ++j )
         {
            Jet jet = l1jet40.at(j);
            for ( int m = 0 ; m < int(l1mu12.size()) ; ++m )
            {
               Muon muon = l1mu12.at(m);
               double deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 )
               {
                  l1muonjets.push_back(jet);
                  l1jetmuons.push_back(muon);
               }
            }
         }
      }
      
      if ( l1jet16.size() >= 1 && l1mu3.size() >= 1 )
      {
         for ( int j = 0 ; j < int(l1jet16.size()) ; ++j )
         {
            Jet jet = l1jet16.at(j);
            for ( int m = 0 ; m < int(l1mu3.size()) ; ++m )
            {
               Muon muon = l1mu3.at(m);
               double deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 )
               {
                  l1muonjetsbtv16.push_back(jet);
                  l1jetmuonsbtv16.push_back(muon);
                  if ( jet.pt() >= 60 )
                  {
                     l1muonjetsbtv60.push_back(jet);
                     l1jetmuonsbtv60.push_back(muon);
                  }
                  if ( jet.pt() >= 120 )
                  {
                     l1muonjetsbtv120.push_back(jet);
                     l1jetmuonsbtv120.push_back(muon);
                  }
               }
            }
         }
      }
      
// --------------------      
      // All hadronic L1 seed fired
      trigger_accept_["L1_DoubjeJet100"]        = ( l1jet100.size() >=2 );
      trigger_accept_["L1_DoubjeJet112"]        = ( l1jet112.size() >=2 );
      trigger_accept_["L1_DoubjeJet100dEta1p6"] = ( l1jet100.size() >=2 && l1jet100deta.size() >= 2 );
      // Semileptonic L1 seed fired
      trigger_accept_["L1_Mu12_DoubleJet40"]        = ( l1mu12.size() >= 1 && l1muonjets.size() >= 1 && l1jet40.size() >= 2 && analysis.triggerResult("HLT_L1SingleMu3_v") ); // the last requirement is to ensure the L1 muon quality selection, not included in the ntuple
      trigger_accept_["L1_Mu12_DoubleJet40dEta1p6"] = ( l1mu12.size() >= 1 && l1muonjets.size() >= 1 && l1jet40.size() >= 2 && l1jet40deta.size() >= 2 && analysis.triggerResult("HLT_L1SingleMu3_v") ); // the last requirement is to ensure the L1 muon quality selection, not included in the ntuple
      trigger_accept_["L1_Mu3_Jet16"]        = ( l1mu3.size() >= 1 && l1muonjetsbtv16.size() >= 1  && l1jet16.size() >= 1 && analysis.triggerResult("HLT_L1SingleMu3_v") ); // the last requirement is to ensure the L1 muon quality selection, not included in the ntuple
      trigger_accept_["L1_Mu3_Jet60"]        = ( l1mu3.size() >= 1 && l1muonjetsbtv60.size() >= 1  && l1jet60.size() >= 1 && analysis.triggerResult("HLT_L1SingleMu3_v") ); // the last requirement is to ensure the L1 muon quality selection, not included in the ntuple
      trigger_accept_["L1_Mu3_Jet120"]       = ( l1mu3.size() >= 1 && l1muonjetsbtv120.size() >= 1 && l1jet120.size() >= 1 && analysis.triggerResult("HLT_L1SingleMu3_v") ); // the last requirement is to ensure the L1 muon quality selection, not included in the ntuple
      trigger_accept_["L1_DoubleJet40"]      = ( l1jet40.size() >=2 );
      
// ======== End of L1 ======== //      
      
// ======== Begin of HLT ======== //  
      
      //
//      bool hltMuJet = false;
//      bool hltPFMuJet = false;
    
// --------------------      
      // Jets collections
      std::shared_ptr< Collection<Jet> >    calojets   = analysis.collection<Jet>("CaloJets");
      std::shared_ptr< Collection<Jet> >    pfjets     = analysis.collection<Jet>("PFJets");
      std::shared_ptr< Collection<JetTag> > jetstags   = analysis.collection<JetTag>("JetsTags");
      // Muons collections 
      std::shared_ptr< Collection<Muon> >   l2muons   = analysis.collection<Muon>("L2Muons");
      std::shared_ptr< Collection<Muon> >   l3muons   = analysis.collection<Muon>("L3Muons");
// --------------------      
// HLT Jet objects
      // Calo
      std::vector<Jet> calojets30;
      std::vector<Jet> calojets100;
      // Calo for b-tagging
      std::vector<Jet> calojets30Selector;
      std::vector<Jet> calojets80Selector;
      // PF
      std::vector<Jet> pfjets40;
      std::vector<Jet> pfjets50;
      std::vector<Jet> pfjets60;
      std::vector<Jet> pfjets70;
      std::vector<Jet> pfjets100;
      std::vector<Jet> pfjets110;
      std::vector<Jet> pfjets112;
      std::vector<Jet> pfjets120;
      std::vector<Jet> pfjets130;
      std::vector<Jet> pfjets140;
      std::vector<Jet> pfjets200;
      std::vector<Jet> pfjets350;
      std::vector<Jet> pfjets40deta;
      std::vector<Jet> pfjets50deta;
      std::vector<Jet> pfjets60deta;
      std::vector<Jet> pfjets70deta;
      std::vector<Jet> pfjets100deta;
      std::vector<Jet> pfjets110deta;
      std::vector<Jet> pfjets120deta;
      std::vector<Jet> pfjets130deta;
// HLT Jet btag      
      std::vector<Jet>  bjets30wp092;
      std::vector<Jet>  bjets100wp092;
      std::vector<Jet>  bjets30wp084;
      std::vector<Jet>  bjets100wp084;
// HLT Muon objects      
      std::vector<Muon> l2muon8;
      std::vector<Muon> l3muon12;
// HLT Muon-Jet objects
      // Calo
      std::vector<Jet>  l2calomuonjets;
      std::vector<Muon> l2calojetmuons;
      std::vector<Jet>  l3calomuonjets;
      std::vector<Muon> l3calojetmuons;
      // PF
      std::vector<Jet>  l2pfmuonjets;
      std::vector<Muon> l2pfjetmuons;
      std::vector<Jet>  l3pfmuonjets40;
      std::vector<Muon> l3pfjetmuons40;
      std::vector<Jet>  l3pfmuonjets50;
      std::vector<Muon> l3pfjetmuons50;
      std::vector<Jet>  l3pfmuonjets60;
      std::vector<Muon> l3pfjetmuons60;
      std::vector<Jet>  l3pfmuonjets70;
      std::vector<Muon> l3pfjetmuons70;
      std::vector<Jet>  l3pfmuonjets100;
      std::vector<Muon> l3pfjetmuons100;
      std::vector<Jet>  l3pfmuonjets200;
      std::vector<Muon> l3pfjetmuons200;
      std::vector<Jet>  l3pfmuonjets350;
      std::vector<Muon> l3pfjetmuons350;
      
// --------------------      
      
      // Jet selection
      for ( int j = 0 ; j < calojets->size() ; ++j )
      {
         Jet jet = calojets->at(j);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= etacut_ ) calojets100.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= etacut_ ) calojets30.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.4 ) calojets30Selector.push_back(jet);
         if ( jet.pt() >= 80  && fabs(jet.eta()) <= 2.4 ) calojets80Selector.push_back(jet);
      }
      
      for ( int j = 0 ; j < pfjets->size() ; ++j )
      {
         Jet jet = pfjets->at(j);
         if ( jet.pt() >= 40  && fabs(jet.eta()) <= etacut_ ) pfjets40.push_back(jet);
         if ( jet.pt() >= 54  && fabs(jet.eta()) <= etacut_ ) pfjets50.push_back(jet);
         if ( jet.pt() >= 62  && fabs(jet.eta()) <= etacut_ ) pfjets60.push_back(jet);
         if ( jet.pt() >= 70  && fabs(jet.eta()) <= etacut_ ) pfjets70.push_back(jet);

                  
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= etacut_ ) pfjets100.push_back(jet);
         if ( jet.pt() >= 112 && fabs(jet.eta()) <= etacut_ ) pfjets112.push_back(jet);
         if ( jet.pt() >= 116 && fabs(jet.eta()) <= etacut_ ) pfjets110.push_back(jet);
         if ( jet.pt() >= 120 && fabs(jet.eta()) <= etacut_ ) pfjets120.push_back(jet);
         if ( jet.pt() >= 130 && fabs(jet.eta()) <= etacut_ ) pfjets130.push_back(jet);
         if ( jet.pt() >= 140 && fabs(jet.eta()) <= etacut_ ) pfjets140.push_back(jet);
         if ( jet.pt() >= 200 && fabs(jet.eta()) <= etacut_ ) pfjets200.push_back(jet);
         if ( jet.pt() >= 350 && fabs(jet.eta()) <= etacut_ ) pfjets350.push_back(jet);
      }
      // L1 delta_eta jets selection
      if ( pfjets40.size() >= 2 )
      {
         for ( int j1 = 0 ; j1 < int(pfjets40.size())-1 ; ++j1 )
         {
            Jet jet1 = pfjets40.at(j1);
            for ( int j2 = j1+1 ; j2 < int(pfjets40.size()) ; ++j2 )
            {
               Jet jet2 = pfjets40.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= detacut_ )
               {
                  // semilep
                  pfjets40deta.push_back(jet1);
                  pfjets40deta.push_back(jet2);
                  if ( jet1.pt() >= 54 && jet2.pt() >= 54 ) // avoid repetition
                  {
                     pfjets50deta.push_back(jet1);
                     pfjets50deta.push_back(jet2);
                  }
                  if ( jet1.pt() >= 62 && jet2.pt() >= 62 ) // avoid repetition
                  {
                     pfjets60deta.push_back(jet1);
                     pfjets60deta.push_back(jet2);
                  }
                  if ( jet1.pt() >= 70 && jet2.pt() >= 70 ) // avoid repetition
                  {
                     pfjets70deta.push_back(jet1);
                     pfjets70deta.push_back(jet2);
                  }
                  // all-hadronic
                  if ( jet1.pt() >= 112 && jet2.pt() >= 112 ) // avoid repetition
                  {
                     pfjets100deta.push_back(jet1);
                     pfjets100deta.push_back(jet2);
                  }
                  if ( jet1.pt() >= 116 && jet2.pt() >= 116 ) // avoid repetition
                  {
                     pfjets110deta.push_back(jet1);
                     pfjets110deta.push_back(jet2);
                  }
                  if ( jet1.pt() >= 120 && jet2.pt() >= 120 ) // avoid repetition
                  {
                     pfjets120deta.push_back(jet1);
                     pfjets120deta.push_back(jet2);
                  }
                  if ( jet1.pt() >= 130 && jet2.pt() >= 130 ) // avoid repetition
                  {
                     pfjets130deta.push_back(jet1);
                     pfjets130deta.push_back(jet2);
                  }
               }
            }
         }
      }
      
      // btagging
      for ( int bj = 0; bj < jetstags->size() ; ++bj )
      {
         JetTag bjet = jetstags->at(bj);
         TLorentzVector bjetp4;
         bjetp4.SetPtEtaPhiM(bjet.pt(),bjet.eta(),bjet.phi(),0.);
         if ( bjet.pt() >= 30 )
         {
            for ( size_t j = 0; j < calojets30.size() ; ++j )
            {
               Jet jet = calojets30.at(j);
               if ( jet.p4().DeltaR(bjetp4) < 0.5 )
               {
                  if ( bjet.btag() >= 0.92 )  bjets30wp092.push_back(jet);
                  if ( bjet.btag() >= 0.84 )  bjets30wp084.push_back(jet);
               }
            }
         }
         if ( bjet.pt() >= 80 && bj < 6 )
         {
            for ( size_t j = 0; j < calojets100.size() ; ++j )
            {
               Jet jet = calojets100.at(j);
               if ( jet.p4().DeltaR(bjetp4) < 0.5 )
               {
                  if ( bjet.btag() >= 0.92 ) bjets100wp092.push_back(jet);
                  if ( bjet.btag() >= 0.84 ) bjets100wp084.push_back(jet);
               }
            }
         }
      }
      
      // Muon selection
      for ( int m = 0 ; m < l2muons->size() ; ++m )
      {
         Muon muon = l2muons->at(m);
         if ( muon.pt() >= 8. && fabs(muon.eta()) <= etacut_ ) l2muon8.push_back(muon);
      }
      
      for ( int m = 0 ; m < l3muons->size() ; ++m )
      {
         Muon muon = l3muons->at(m);
         if ( muon.pt() >= 12. && fabs(muon.eta()) <= etacut_ ) l3muon12.push_back(muon);
      }
      
      // Mu-Jet selection
      // Calo BTagMu L2.5: Calo Jets + L2 Muons
      if ( calojets30.size() >= 1 && l2muon8.size() >= 1 )
      {
         for ( int j = 0; j < int(calojets30.size()) ; ++j )
         {
            Jet jet = calojets30.at(j);
            for ( int m = 0; m < int(l2muon8.size()); ++m )
            {
               Muon muon = l2muon8.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.5 ) 
               {   
                  l2calomuonjets.push_back(jet);
                  l2calojetmuons.push_back(muon);
               }
            }
         }
      }
      // Calo BTagMu L3: Calo Jets + L3 Muons
      if ( calojets30.size() >= 1 && l3muon12.size() >= 1 )
      {
         for ( int j = 0; j < int(calojets30.size()) ; ++j )
         {
            Jet jet = calojets30.at(j);
            for ( int m = 0; m < int(l3muon12.size()); ++m )
            {
               Muon muon = l3muon12.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.5 ) 
               {   
                  l3calomuonjets.push_back(jet);
                  l3calojetmuons.push_back(muon);
               }
            }
         }
      }
      
      // PF BTagMu L2.5: PF Jets + L2 Muons
      if ( pfjets40.size() >= 1 && l2muon8.size() >= 1 )
      {
         for ( int j = 0; j < int(pfjets40.size()) ; ++j )
         {
            Jet jet = pfjets40.at(j);
            for ( int m = 0; m < int(l2muon8.size()); ++m )
            {
               Muon muon = l2muon8.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 ) 
               {   
                  l2pfmuonjets.push_back(jet);
                  l2pfjetmuons.push_back(muon);
               }
            }
         }
      }
      // PF BTagMu L3: PF Jets + L3 Muons
      if ( pfjets40.size() >= 1 && l3muon12.size() >= 1 )
      {
         for ( int j = 0; j < int(pfjets40.size()) ; ++j )
         {
            Jet jet = pfjets40.at(j);
            for ( int m = 0; m < int(l3muon12.size()); ++m )
            {
               Muon muon = l3muon12.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 ) 
               {   
                  l3pfmuonjets40.push_back(jet);
                  l3pfjetmuons40.push_back(muon);
                  if ( jet.pt() >= 54 )
                  {
                     l3pfmuonjets50.push_back(jet);
                     l3pfjetmuons50.push_back(muon);
                  }
                  if ( jet.pt() >= 62 )
                  {
                     l3pfmuonjets60.push_back(jet);
                     l3pfjetmuons60.push_back(muon);
                  }
                  if ( jet.pt() >= 70 )
                  {
                     l3pfmuonjets70.push_back(jet);
                     l3pfjetmuons70.push_back(muon);
                  }
                  if ( jet.pt() >= 100 )
                  {
                     l3pfmuonjets100.push_back(jet);
                     l3pfjetmuons100.push_back(muon);
                  }
                  if ( jet.pt() >= 200 )
                  {
                     l3pfmuonjets200.push_back(jet);
                     l3pfjetmuons200.push_back(muon);
                  }
                  if ( jet.pt() >= 350 )
                  {
                     l3pfmuonjets350.push_back(jet);
                     l3pfjetmuons350.push_back(muon);
                  }
               }
            }
         }
      }
//      hltMuJet   = ( l2calomuonjets.size() >= 1 && l3calomuonjets.size() >= 1 && l2pfmuonjets.size() >= 1 && l3pfmuonjets.size() >= 1 );
//      hltPFMuJet = ( l2pfmuonjets.size() >= 1 && l3pfmuonjets.size() >= 1 );

      // All hadronic
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"]  = ( calojets100.size() >= 2 && bjets100wp084.size() >= 2 && pfjets100.size() >= 2 && pfjets100deta.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"]  = ( calojets100.size() >= 2 && bjets100wp092.size() >= 2 && pfjets112.size() >= 2 && pfjets100deta.size() >= 2 );
      // monitor
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 1 && pfjets100.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 1 && pfjets140.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 1 && pfjets200.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 1 && pfjets350.size() >= 2 );
                   
      // backup
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 2 && pfjets110.size() >= 2 && pfjets110deta.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 2 && pfjets120.size() >= 2 && pfjets120deta.size() >= 2 );
      trigger_accept_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130"] = ( calojets100.size() >= 2 && bjets100wp092.size() >= 2 && pfjets130.size() >= 2 && pfjets130deta.size() >= 2 );
                                                               
      
      // Semileptonic
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] = ( calojets30.size() >= 2 && bjets30wp084.size() >= 2 && pfjets40.size() >= 2 && pfjets40deta.size() >= 2 && l3pfmuonjets40.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 2 && pfjets40.size() >= 2 && pfjets40deta.size() >= 2 && l3pfmuonjets40.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      // monitor
      trigger_accept_["HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40"]         = ( calojets30.size() >= 2 && bjets30wp092.size() >= 1 && pfjets40.size()  >= 2 );
      trigger_accept_["HLT_1CaloJets30_Muon12_1PFJets40"]                  = ( calojets30.size() >= 1 && pfjets40.size()  >= 1 && l3pfmuonjets40.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40"]  = ( calojets30.size() >= 2 && bjets30wp092.size() >= 1 && pfjets40.size()  >= 2 && l3pfmuonjets40.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 1 && pfjets100.size() >= 2 && l3pfmuonjets100.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 1 && pfjets200.size() >= 2 && l3pfmuonjets200.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 1 && pfjets350.size() >= 2 && l3pfmuonjets350.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      // backup
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 2 && pfjets50.size() >= 2 && pfjets50deta.size() >= 2 && l3pfmuonjets50.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 2 && pfjets60.size() >= 2 && pfjets60deta.size() >= 2 && l3pfmuonjets60.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );
      trigger_accept_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70"] = ( calojets30.size() >= 2 && bjets30wp092.size() >= 2 && pfjets70.size() >= 2 && pfjets70deta.size() >= 2 && l3pfmuonjets70.size() >= 1 && l2muon8.size() >= 1 && l3muon12.size() >= 1 );    
      
      
// --------------------      
      // HLT Path - All had
      if ( trigger_accept_["L1_DoubjeJet112"] )
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
      if ( trigger_accept_["L1_Mu3_Jet16"] && trigger_accept_["HLT_1CaloJets30_Muon12_1PFJets40"] )
      {
         h1_["HLT_1CaloJets30_Muon12_1PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_1CaloJets30_Muon12_1PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
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
