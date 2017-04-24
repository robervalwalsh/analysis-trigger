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
   
   int ncoll = 2200;
   int ncollFix = -1;

   xsections_ = CrossSections();
   
   triggers_.push_back("ZeroBias");
   
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40");
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV088_2PFJets40");
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40");
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV094_2PFJets40");
   
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40");                  // v2
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1");           // v3
   triggers_.push_back("HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5");   // v4
   
   triggers_.push_back("HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100");
   triggers_.push_back("HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100");
   triggers_.push_back("HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1");
   triggers_.push_back("HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5");

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
   if ( ! isData_ )  // temporary
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
      if ( ! isData_ ) //temporary
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
      
      // L1 Jets
      std::shared_ptr< Collection<Jet> > l1jets = analysis.collection<Jet>("L1Jets");
      std::vector<Jet> l1jet32eta2p3;
      std::vector<Jet> l1jet100eta2p3;
      for ( int j = 0; j < l1jets->size() ; ++j )
      {
         Jet l1j = l1jets->at(j);
         if ( l1j.pt() >= 32  && fabs(l1j.eta()) <= 2.3 ) l1jet32eta2p3.push_back(l1j);
         if ( l1j.pt() >= 100 && fabs(l1j.eta()) <= 2.3 ) l1jet100eta2p3.push_back(l1j);
      }
      
      // L1 Muons
      std::shared_ptr< Collection<Muon> > l1muons = analysis.collection<Muon>("L1Muons");
      std::vector<Muon> l1mu10eta2p3;
      for ( int m = 0; m < l1muons->size() ; ++m )
      {
         Muon l1m = l1muons->at(m);
         if ( l1m.pt() >= 10  && fabs(l1m.eta()) <= 2.3 ) l1mu10eta2p3.push_back(l1m);
      }
      
      
      // L1 All hadronic seed
      bool l1Allhad = false;
      if ( l1jet100eta2p3.size() >= 2 )
      {
         for ( int j1 = 0 ; j1 < int(l1jet100eta2p3.size())-1 ; ++j1 )
         {
            Jet jet1 = l1jet100eta2p3.at(j1);
            for ( int j2 = j1+1 ; j2 < int(l1jet100eta2p3.size()) ; ++j2 )
            {
               Jet jet2 = l1jet100eta2p3.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) l1Allhad = true;
            }
         }
      }
      
      // L1 semileptonic seed
      bool l1Semilep = false;
      // L1 jet-mu
      bool l1JetMu = false;
      if ( l1jet32eta2p3.size() >= 1 && l1mu10eta2p3.size() >= 1 )
      {
         for ( int j = 0 ; j < int(l1jet32eta2p3.size()) ; ++j )
         {
            Jet jet = l1jet32eta2p3.at(j);
            for ( int m = 0 ; m < int(l1mu10eta2p3.size()) ; ++m )
            {
               Muon muon = l1mu10eta2p3.at(m);
               double deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 ) l1JetMu = true;
            }
         }
      }
      
      if ( l1jet32eta2p3.size() >= 2 && l1JetMu )
      {
         for ( int j1 = 0 ; j1 < int(l1jet32eta2p3.size())-1 ; ++j1 )
         {
            Jet jet1 = l1jet32eta2p3.at(j1);
            for ( int j2 = j1+1 ; j2 < int(l1jet32eta2p3.size()) ; ++j2 )
            {
               Jet jet2 = l1jet32eta2p3.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) l1Semilep = true;
            }
         }
      }
      // HLT
      // Jets   
      std::shared_ptr< Collection<Jet> >    calojets   = analysis.collection<Jet>("CaloJets");
      std::shared_ptr< Collection<Jet> >    pfjets     = analysis.collection<Jet>("PFJets");
      std::shared_ptr< Collection<JetTag> > jetstags   = analysis.collection<JetTag>("JetsTags");
      // Muons   
      std::shared_ptr< Collection<Muon> >   l2muons   = analysis.collection<Muon>("L2Muons");
      std::shared_ptr< Collection<Muon> >   l3muons   = analysis.collection<Muon>("L3Muons");

      
      // Jet selection
      std::vector<Jet>  calojets30eta2p3;
      std::vector<Jet>  calojets30eta2p1;
      std::vector<Jet>  calojets30Selector;
      std::vector<Jet>  calojets30eta2p1Selector;
      std::vector<Jet>  calojets100eta2p3;
      std::vector<Jet>  calojets100eta2p1;
      std::vector<Jet>  calojets80Selector;
      std::vector<Jet>  calojets80eta2p1Selector;
      for ( int j = 0 ; j < calojets->size() ; ++j )
      {
         Jet jet = calojets->at(j);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= 2.3 ) calojets100eta2p3.push_back(jet);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= 2.1 ) calojets100eta2p1.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.3 ) calojets30eta2p3.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.1 ) calojets30eta2p1.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.4 ) calojets30Selector.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.1 ) calojets30eta2p1Selector.push_back(jet);
         if ( jet.pt() >= 80  && fabs(jet.eta()) <= 2.4 ) calojets80Selector.push_back(jet);
         if ( jet.pt() >= 80  && fabs(jet.eta()) <= 2.1 ) calojets80eta2p1Selector.push_back(jet);
      }
      std::vector<Jet> pfjets100eta2p3;
      std::vector<Jet> pfjets100eta2p1;
      std::vector<Jet> pfjets40eta2p3;
      std::vector<Jet> pfjets50eta2p3;
      std::vector<Jet> pfjets40eta2p1;
      std::vector<Jet> pfjets50eta2p1;
      
      for ( int j = 0 ; j < pfjets->size() ; ++j )
      {
         Jet jet = pfjets->at(j);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= 2.3 ) pfjets100eta2p3.push_back(jet);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= 2.1 ) pfjets100eta2p1.push_back(jet);
         if ( jet.pt() >= 40  && fabs(jet.eta()) <= 2.3 ) pfjets40eta2p3.push_back(jet);
         if ( jet.pt() >= 50  && fabs(jet.eta()) <= 2.3 ) pfjets50eta2p3.push_back(jet);
         if ( jet.pt() >= 40  && fabs(jet.eta()) <= 2.1 ) pfjets40eta2p1.push_back(jet);
         if ( jet.pt() >= 50  && fabs(jet.eta()) <= 2.1 ) pfjets50eta2p1.push_back(jet);
      }
      
      // Muon selection
      std::vector<Muon> l2muon8eta2p3;
      std::vector<Muon> l3muon12eta2p3;
      std::vector<Muon> l2muon8eta2p1;
      std::vector<Muon> l3muon12eta2p1;
      for ( int m = 0 ; m < l2muons->size() ; ++m )
      {
         Muon muon = l2muons->at(m);
         if ( muon.pt() >= 8. && fabs(muon.eta()) <= 2.3 ) l2muon8eta2p3.push_back(muon);
         if ( muon.pt() >= 8. && fabs(muon.eta()) <= 2.1 ) l2muon8eta2p1.push_back(muon);
      }
      
      for ( int m = 0 ; m < l3muons->size() ; ++m )
      {
         Muon muon = l3muons->at(m);
         if ( muon.pt() >= 12. && fabs(muon.eta()) <= 2.3 ) l3muon12eta2p3.push_back(muon);
         if ( muon.pt() >= 12. && fabs(muon.eta()) <= 2.1 ) l3muon12eta2p1.push_back(muon);
      }
      
      // btagging
      std::vector<JetTag>  jetstags30wp084;
      std::vector<JetTag>  jetstags30wp088;
      std::vector<JetTag>  jetstags30wp092;
      std::vector<JetTag>  jetstags30wp094;
      std::vector<JetTag>  jetstags30eta2p1wp092;
      
      std::vector<JetTag>  jetstags100wp084;
      std::vector<JetTag>  jetstags100wp092;
      std::vector<JetTag>  jetstags100eta2p1wp092;

      for ( int jt = 0; jt < jetstags->size() ; ++jt )
      {
         JetTag jet = jetstags->at(jt);
         if ( jet.pt() >= 30 )
         {
            if ( jet.btag() >= 0.84 ) jetstags30wp084.push_back(jet);
            if ( jet.btag() >= 0.88 ) jetstags30wp088.push_back(jet);
            if ( jet.btag() >= 0.92 ) jetstags30wp092.push_back(jet);
            if ( jet.btag() >= 0.94 ) jetstags30wp094.push_back(jet);
            if ( jet.btag() >= 0.92 && fabs(jet.btag()) <= 2.1 ) jetstags30eta2p1wp092.push_back(jet);
         }
         if ( jet.pt() >= 80 )
         {
            if ( jet.btag() >= 0.84 ) jetstags100wp084.push_back(jet);
            if ( jet.btag() >= 0.92 ) jetstags100wp092.push_back(jet);
            if ( jet.btag() >= 0.92 && fabs(jet.btag()) <= 2.1 ) jetstags100eta2p1wp092.push_back(jet);
         }
      }
      
      // Mu-Jet
      bool hltMuJet = false;
      std::vector<Jet>  muonjets40eta2p3;
      std::vector<Muon> jetmuons40eta2p3;
      if ( l1Semilep && pfjets40eta2p3.size() >= 1 && l3muon12eta2p3.size() >= 1 )
      {
         for ( int j = 0; j < int(pfjets40eta2p3.size()) ; ++j )
         {
            Jet jet = pfjets40eta2p3.at(j);
            for ( int m = 0; m < int(l3muon12eta2p3.size()); ++m )
            {
               Muon muon = l3muon12eta2p3.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 ) 
               {   
                  muonjets40eta2p3.push_back(jet);
                  jetmuons40eta2p3.push_back(muon);
               }
            }
         }
      }
      hltMuJet = ( muonjets40eta2p3.size() >= 1 );
      
      bool hltMuJeteta2p1 = false;
      std::vector<Jet>  muonjets40eta2p1;
      std::vector<Muon> jetmuons40eta2p1;
      if ( l1Semilep && pfjets40eta2p1.size() >= 1 && l3muon12eta2p1.size() >= 1 )
      {
         for ( int j = 0; j < int(pfjets40eta2p1.size()) ; ++j )
         {
            Jet jet = pfjets40eta2p1.at(j);
            for ( int m = 0; m < int(l3muon12eta2p1.size()); ++m )
            {
               Muon muon = l3muon12eta2p1.at(m);
               float deltaR = jet.p4().DeltaR(muon.p4());
               if ( deltaR <= 0.4 ) 
               {   
                  muonjets40eta2p1.push_back(jet);
                  jetmuons40eta2p1.push_back(muon);
               }
            }
         }
      }
      hltMuJeteta2p1 = ( muonjets40eta2p1.size() >= 1 );
      
      // HLT Path - Semileptonic
      bool hltSemilep084 = false;
      bool hltSemilep088 = false;
      bool hltSemilep092 = false;
      bool hltSemilep094 = false;
      bool hltSemilep092v2 = false;
      bool hltSemilep092v3 = false;
      bool hltSemilep092v4 = false;
      if ( l1Semilep && calojets30eta2p3.size() >= 2 && pfjets40eta2p3.size() >= 2 && l3muon12eta2p3.size() >= 1 && hltMuJet )
      {
         // delta eta
         bool deltaEta1p6 = false;
         for ( int j1 = 0; j1 < int(pfjets40eta2p3.size()-1); ++j1 )
         {
            Jet jet1 = pfjets40eta2p3.at(j1);
            for ( int j2 = j1+1; j2 < int(pfjets40eta2p3.size()); ++j2 )
            {
               Jet jet2 = pfjets40eta2p3.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) deltaEta1p6 = true;
            }
         }
         if ( jetstags30wp084.size() >= 2 && deltaEta1p6 )    hltSemilep084 = true;
         if ( jetstags30wp088.size() >= 2 && deltaEta1p6 )    hltSemilep088 = true;
         if ( jetstags30wp092.size() >= 2 && deltaEta1p6 )    hltSemilep092 = true;
         if ( jetstags30wp094.size() >= 2 && deltaEta1p6 )    hltSemilep094 = true;
         
         if ( jetstags30wp092.size() >= 2 && deltaEta1p6 && pfjets50eta2p3.size() >= 1 )    hltSemilep092v2 = true;
      }
      
      if ( l1Semilep && calojets30eta2p1.size() >= 2 && pfjets40eta2p1.size() >= 2 && l3muon12eta2p1.size() >= 1 && hltMuJeteta2p1 )
      {
         // delta eta
         bool deltaEta1p6 = false;
         bool deltaEta1p5 = false;
         for ( int j1 = 0; j1 < int(pfjets40eta2p1.size()-1); ++j1 )
         {
            Jet jet1 = pfjets40eta2p1.at(j1);
            for ( int j2 = j1+1; j2 < int(pfjets40eta2p1.size()); ++j2 )
            {
               Jet jet2 = pfjets40eta2p1.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) deltaEta1p6 = true;
               if ( deltaEta <= 1.5 ) deltaEta1p5 = true;
            }
         }
         if ( jetstags30eta2p1wp092.size() >= 2 && deltaEta1p6 && pfjets50eta2p1.size() >= 1 )    hltSemilep092v3 = true;
         if ( jetstags30eta2p1wp092.size() >= 2 && deltaEta1p5 && pfjets50eta2p1.size() >= 1 )    hltSemilep092v4 = true;
      }
      
      if ( hltSemilep084 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltSemilep088 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV088_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV088_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltSemilep092 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltSemilep094 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV094_2PFJets40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV094_2PFJets40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      
      if ( hltSemilep092v2 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_psw"] -> Fill(nPileup, prescale/ncoll);   
      }

      if ( hltSemilep092v3 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      
      if ( hltSemilep092v4 )
      {
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      
      // HLT Path - All had
      bool hltAllhad084 = false;
      bool hltAllhad092 = false;
      bool hltAllhad092v3 = false;
      bool hltAllhad092v4 = false;
      if ( l1Allhad && calojets100eta2p3.size() >= 2 && pfjets100eta2p3.size() >= 2 )
      {
         // delta eta
         bool deltaEta1p6 = false;
         for ( int j1 = 0; j1 < int(pfjets100eta2p3.size()-1); ++j1 )
         {
            Jet jet1 = pfjets100eta2p3.at(j1);
            for ( int j2 = j1+1; j2 < int(pfjets100eta2p3.size()); ++j2 )
            {
               Jet jet2 = pfjets100eta2p3.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) deltaEta1p6 = true;
            }
         }
         
         if ( jetstags100wp084.size() >= 2 && deltaEta1p6 )    hltAllhad084 = true;
         if ( jetstags100wp092.size() >= 2 && deltaEta1p6 )    hltAllhad092 = true;
      }
      // $$$$$$$$$$$$$$
      if ( l1Allhad && calojets100eta2p1.size() >= 2 && pfjets100eta2p1.size() >= 2  )
      {
         // delta eta
         bool deltaEta1p6 = false;
         bool deltaEta1p5 = false;
         for ( int j1 = 0; j1 < int(pfjets100eta2p1.size()-1); ++j1 )
         {
            Jet jet1 = pfjets100eta2p1.at(j1);
            for ( int j2 = j1+1; j2 < int(pfjets100eta2p1.size()); ++j2 )
            {
               Jet jet2 = pfjets100eta2p1.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= 1.6 ) deltaEta1p6 = true;
               if ( deltaEta <= 1.5 ) deltaEta1p5 = true;
            }
         }
         if ( jetstags100eta2p1wp092.size() >= 2 && deltaEta1p6 )    hltAllhad092v3 = true;
         if ( jetstags100eta2p1wp092.size() >= 2 && deltaEta1p5 )    hltAllhad092v4 = true;
      }
      
      
      if ( hltAllhad084 )
      {
         h1_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltAllhad092 )
      {
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltAllhad092v3 )
      {
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
      if ( hltAllhad092v4 )
      {
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5"] -> Fill(nPileup, 1./ncoll);   
         h1_["HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5_psw"] -> Fill(nPileup, prescale/ncoll);   
      }
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
