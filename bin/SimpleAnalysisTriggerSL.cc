#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Core/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::core;


// =============================================================================================   
int main(int argc, char * argv[])
{
   bool isMC = false;
   bool isbbb = false;
   std::string inputList = "rootFilelist_v2.txt";
   std::string outputRoot = "histograms_ntuples_purity_bbnb.root";
   std::string json = "json_2016H_v2.txt";
   int bWP = 2;
   float btagcut[3] = {0.46,0.8434,0.935};
   // Cuts                                         // <<<<===== CMSDAS
   float ptmin[3]   = { 40.0, 40.0, 30.0 };
   float muptmin    = 12;
   float etamax[3]  = {   2.2,   2.2 , 2.2 };
   float muetamax   = 2.2;
   float btagmin[3] = { btagcut[bWP], btagcut[bWP], btagcut[bWP]};
   float nonbtag    = 0.46;
   float dRmin      = 1.;
   float detamax    = 1.55;
   float etacut_ = 2.3;
   float detacut_ = 1.6;
   
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   Analysis analysis(inputList,"MssmHbbTrigger/Events/EventInfo");
   
   // Jets
   analysis.addTree<Jet> ("Jets","MssmHbbTrigger/Events/slimmedJetsPuppi");
   analysis.addTree<Jet> ("L1Jets","MssmHbbTrigger/Events/l1tJets");
   analysis.addTree<Jet> ("CaloJets","MssmHbbTrigger/Events/hltAK4CaloJetsCorrectedIDPassed");
   analysis.addTree<Jet> ("PFJets","MssmHbbTrigger/Events/hltAK4PFJetsCorrected");
   analysis.addTree<Jet> ("hltJetsIdLoose","MssmHbbTrigger/Events/hltAK4PFJetsLooseIDCorrected");
   
   // Muons
   analysis.addTree<Muon> ("L1Muons","MssmHbbTrigger/Events/l1tMuons");
   analysis.addTree<Muon> ("L2Muons","MssmHbbTrigger/Events/hltL2MuonCandidates");
   analysis.addTree<Muon> ("L3Muons","MssmHbbTrigger/Events/hltL3MuonCandidates");
   analysis.addTree<Muon> ("Muons","MssmHbbTrigger/Events/slimmedMuons");
   
   // Vertices
   analysis.addTree<Vertex> ("FastPV","MssmHbbTrigger/Events/hltFastPrimaryVertex");
   analysis.addTree<Vertex> ("FastPVPixel","MssmHbbTrigger/Events/hltFastPVPixelVertices");
   analysis.addTree<Vertex> ("offlinePV","MssmHbbTrigger/Events/offlineSlimmedPrimaryVertices");
   
   // BTag
   analysis.addTree<JetTag> ("JetsTags","MssmHbbTrigger/Events/hltCombinedSecondaryVertexBJetTagsCalo");

   
   std::vector<std::string> triggerObjects;
   triggerObjects.push_back("hltL1sDoubleJetC100");
   triggerObjects.push_back("hltDoubleJetsC100");
   triggerObjects.push_back("hltBTagCaloCSVp014DoubleWithMatching");
   triggerObjects.push_back("hltDoublePFJetsC100");
   triggerObjects.push_back("hltDoublePFJetsC100MaxDeta1p6");

   for ( auto & obj : triggerObjects )
      analysis.addTree<TriggerObject> (obj,Form("MssmHbbTrigger/Events/hltTriggerSummaryAOD/%s",obj.c_str()));
   
   analysis.triggerResults("MssmHbbTrigger/Events/TriggerResults");
   std::string hltPath = "HLT_BTagMu_DiJet20_Mu5_v";
   
   std::cout << "oioi" << std::endl;
   
   if( !isMC ) analysis.processJsonFile(json);
   
   TFile hout(outputRoot.c_str(),"recreate");
   
   std::map<std::string, TH1F*> h1;
   h1["n"]        = new TH1F("n" , "" , 30, 0, 30);
   h1["n_csv"]    = new TH1F("n_csv" , "" , 30, 0, 30);
   h1["n_ptmin20"]= new TH1F("n_ptmin20" , "" , 30, 0, 30);
   h1["n_ptmin20_csv"] = new TH1F("n_ptmin20_csv" , "" , 30, 0, 30);
   for ( int i = 0 ; i < 3 ; ++i )
   {
      h1[Form("pt_%i",i)]         = new TH1F(Form("pt_%i",i) , "" , 100, 0, 1000);
      h1[Form("eta_%i",i)]        = new TH1F(Form("eta_%i",i) , "" , 100, -5, 5);
      h1[Form("phi_%i",i)]        = new TH1F(Form("phi_%i",i) , "" , 100, -4, 4);
      h1[Form("btag_%i",i)]       = new TH1F(Form("btag_%i",i) , "" , 100, 0, 1);
      
      h1[Form("pt_%i_csv",i)]     = new TH1F(Form("pt_%i_csv",i) , "" , 100, 0, 1000);
      h1[Form("eta_%i_csv",i)]    = new TH1F(Form("eta_%i_csv",i) , "" , 100, -5, 5);
      h1[Form("phi_%i_csv",i)]    = new TH1F(Form("phi_%i_csv",i) , "" , 100, -4, 4);
      h1[Form("btag_%i_csv",i)]   = new TH1F(Form("btag_%i_csv",i) , "" , 100, 0, 1);
   }
   h1["m12"]     = new TH1F("m12"     , "" , 50, 0, 1000);
   h1["m12_csv"] = new TH1F("m12_csv" , "" , 50, 0, 1000);
   
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: kinematics
   // 3: matched to online
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
   int nsel[10] = { };
   int nmatch[10] = { };
   
   std::string prevFile = "";
//   int mycounter = 0;
   for ( int i = 0 ; i < analysis.size() ; ++i )
//   for ( int i = 0 ; i < 500000 ; ++i )
   {
      int njets = 0;
      int njets_csv = 0;
      bool goodEvent = true;
      
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      if (! isMC )
      {
//          if ( analysis.run() == 279975)
//          {
//             if ( prevFile == "" || prevFile != analysis.fileName() )
//             {
//                prevFile = analysis.fileName();
//                std::cout << prevFile << std::endl;
//             }
//          }
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      // primary vertex
      std::shared_ptr< Collection<Vertex> > fastPVs = analysis.collection<Vertex>("FastPV");
      std::shared_ptr< Collection<Vertex> > fastPVsPixel = analysis.collection<Vertex>("FastPVPixel");
      std::shared_ptr< Collection<Vertex> > offPVs = analysis.collection<Vertex>("offlinePV");
      if ( fastPVs -> size() < 1 || fastPVsPixel -> size() < 1 ) continue;
      Vertex fastpv = fastPVs->at(0);
      if ( !( ! fastpv.fake() && fastpv.ndof() > 0 && fabs(fastpv.z()) <= 25 && fastpv.rho() <= 2 ) ) continue;
      Vertex fastpvpix = fastPVsPixel->at(0);
      if ( !( ! fastpvpix.fake() && fastpvpix.ndof() > 0 && fabs(fastpvpix.z()) <= 25 && fastpvpix.rho() <= 2 ) ) continue;
      
      int triggerFired = analysis.triggerResult(hltPath);
      if ( !triggerFired ) continue;
      
      // L1
// --------------------      
// L1 Jets collection
      std::shared_ptr< Collection<Jet> > l1jets = analysis.collection<Jet>("L1Jets");
// L1 Muons collection
      std::shared_ptr< Collection<Muon> > l1muons = analysis.collection<Muon>("L1Muons");
// L1 jets objects
      std::vector<Jet> l1jet40;
// L1 delta_eta jets objects
      std::vector<Jet> l1jet40deta;
// L1 Muons objects
      std::vector<Muon> l1mu12;
// L1 jet-mu objects
      std::vector<Jet> l1muonjets;
      std::vector<Muon> l1jetmuons;
// --------------------      
      // L1 jets selection
      for ( int j = 0; j < l1jets->size() ; ++j )
      {
         Jet l1j = l1jets->at(j);
         if ( l1j.pt() >= 40  && fabs(l1j.eta()) <= etacut_ ) l1jet40.push_back(l1j);
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
               }
            }
         }
      }
      // L1 muon selection
      for ( int m = 0; m < l1muons->size() ; ++m )
      {
         Muon l1m = l1muons->at(m);
         if ( l1m.pt() >= 12  && fabs(l1m.eta()) <= etacut_ ) l1mu12.push_back(l1m);
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
      
      if ( ! ( l1mu12.size() >= 1 && l1muonjets.size() >= 1 && l1jet40.size() >= 2 && l1jet40deta.size() >= 2 && analysis.triggerResult("HLT_L1SingleMu3_v") ) ) continue;
      
      // HLT
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
      // Calo for b-tagging
      std::vector<Jet> calojets30Selector;
      // PF
      std::vector<Jet> pfjets40;
      std::vector<Jet> pfjets40deta;
// HLT Jet btag      
      std::vector<JetTag>  bjets30wp092;
      std::vector<JetTag>  bjets30wp084;
// HLT Muon objects      
      std::vector<Muon> l2muon8;
      std::vector<Muon> l3muon12;
// HLT Muon-Jet objects
      // PF
      std::vector<Jet>  l3pfmuonjets40;
      std::vector<Muon> l3pfjetmuons40;
      
      // -------------      
            
      // Jet selection
      for ( int j = 0 ; j < calojets->size() ; ++j )
      {
         Jet jet = calojets->at(j);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= etacut_ ) calojets30.push_back(jet);
         if ( jet.pt() >= 30  && fabs(jet.eta()) <= 2.4 ) calojets30Selector.push_back(jet);
      }
      if ( calojets30.size() < 2 ) continue;
      
      for ( int j = 0 ; j < pfjets->size() ; ++j )
      {
         Jet jet = pfjets->at(j);
         if ( jet.pt() >= 40  && fabs(jet.eta()) <= etacut_ ) pfjets40.push_back(jet);
      }
      if ( pfjets40.size() < 2 ) continue;
      
      
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
               }
            }
         }
      }
      
      if ( pfjets40deta.size() < 2 ) continue;
      
      // Muon selection
      for ( int m = 0 ; m < l2muons->size() ; ++m )
      {
         Muon muon = l2muons->at(m);
         if ( muon.pt() >= 8. && fabs(muon.eta()) <= etacut_ ) l2muon8.push_back(muon);
      }
      
      if ( l2muon8.size() < 1 ) continue;
      
      for ( int m = 0 ; m < l3muons->size() ; ++m )
      {
         Muon muon = l3muons->at(m);
         if ( muon.pt() >= 12. && fabs(muon.eta()) <= etacut_ ) l3muon12.push_back(muon);
      }
      
      if ( l3muon12.size() < 1 ) continue;
      
      // Mu-Jet selection
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
               }
            }
         }
      }
      
      if ( l3pfmuonjets40.size() < 1 ) continue;
      
      // online btagging
      for ( int bj = 0; bj < jetstags->size() ; ++bj )
      {
         JetTag bjet = jetstags->at(bj);
         if (bjet.pt() >= 30 && bj < 8 )
         {
            if ( bjet.btag() >= 0.92 ) bjets30wp092.push_back(bjet);
         }
      }
      if ( bjets30wp092.size() < 2 ) continue;
      
//      ++mycounter;
//      std::cout  << mycounter << std::endl;
      
      ++nsel[0];
      
      // match offline to online
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects,0.5);
      Vertex pv = offPVs->at(0);
      if ( !( ! pv.fake() && pv.ndof() > 4 && fabs(pv.z()) <= 24 && pv.rho() <= 2 ) ) continue;
      
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         if ( slimmedJets->at(j).idLoose() ) selectedJets.push_back(&slimmedJets->at(j));
      }
      if ( selectedJets.size() < 3 ) continue;
      
      ++nsel[1];
      
      // Kinematic selection - 3 leading jets
      for ( int j = 0; j < 3; ++j )
      {
         Jet * jet = selectedJets[j];
         if ( jet->pt() < ptmin[j] || fabs(jet->eta()) > etamax[j] )
         {
            goodEvent = false;
            break;
         }
      }
      
      if ( ! goodEvent ) continue;

//******      PAREI AQUI TEM QUE COLOCAR UM MUON NUM DOS 2 JETS
      // Muon - std::shared_ptr< Collection<Jet> >
      auto slimmedMuons = analysis.collection<Muon>("Muons");
      goodEvent = false;
      for ( int j = 0 ; j < slimmedMuons->size() ; ++j )
      {
         Muon muon = slimmedMuons->at(j);
         if ( muon.pt() > muptmin && fabs(muon.eta()) < muetamax )
         {
            if ( muon.p4().DeltaR(selectedJets[0]->p4()) < 0.4 || muon.p4().DeltaR(selectedJets[1]->p4()) < 0.4 )
            {
               goodEvent = true;
               break;
            }
         }
      
      }

      if ( ! goodEvent ) continue;
      
            
      ++nsel[2];
      
      for ( int j1 = 0; j1 < 2; ++j1 )
      {
         const Jet & jet1 = *selectedJets[j1];
         for ( int j2 = j1+1; j2 < 3; ++j2 )
         {
            const Jet & jet2 = *selectedJets[j2];
            if ( jet1.deltaR(jet2) < dRmin ) goodEvent = false;
         }
      }
      
      if ( ! goodEvent ) continue;
      
      ++nsel[3];
      
      if ( fabs(selectedJets[0]->eta() - selectedJets[1]->eta()) > detamax ) continue;
      
      ++nsel[4];
      
      
      // Fill histograms of kinematic passed events
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if ( selectedJets[j]->pt() < 20. ) continue;
         ++njets;
      }
      
      h1["n"] -> Fill(selectedJets.size());
      h1["n_ptmin20"] -> Fill(njets);
      for ( int j = 0; j < 3; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i",j)]  -> Fill(jet->phi());
         h1[Form("btag_%i",j)] -> Fill(jet->btag());
         
         if ( j < 2 && jet->btag() < btagmin[j] )     goodEvent = false;
         if ( ! isbbb )
         {
            if ( j == 2 && jet->btag() > nonbtag )    goodEvent = false; 
         }
         else
         {
            if ( j == 2 && jet->btag() < btagmin[j] ) goodEvent = false; 
         }
      }
      
      h1["m12"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
      
      if ( ! goodEvent ) continue;
      
      ++nsel[5];
      
      // Is matched?
//       bool matched[10] = {true,true,true,true,true,true,true,true,true,true};
//       for ( int j = 0; j < 2; ++j )
//       {
//          Jet * jet = selectedJets[j];
// //         for ( auto & obj : triggerObjects )   matched = (matched && jet->matched(obj));
//          for ( size_t io = 0; io < triggerObjects.size() ; ++io )
//          {       
//             if ( ! jet->matched(triggerObjects[io]) ) matched[io] = false;
//          }
//       }
//       
//       for ( size_t io = 0; io < triggerObjects.size() ; ++io )
//       {
//          if ( matched[io] ) ++nmatch[io];
//          goodEvent = ( goodEvent && matched[io] );
//       }
//       
//       if ( ! goodEvent ) continue;
      
      ++nsel[6];
     
      // Fill histograms of passed bbnb btagging selection
      for ( int j = 0 ; j < (int)selectedJets.size() ; ++j )
      {
         if ( selectedJets[j]->pt() < 20. ) continue;
         ++njets_csv;
      }
      h1["n_csv"] -> Fill(selectedJets.size());
      h1["n_ptmin20_csv"] -> Fill(njets_csv);
      for ( int j = 0; j < 3; ++j )
      {
         Jet * jet = selectedJets[j];
         h1[Form("pt_%i_csv",j)]   -> Fill(jet->pt());
         h1[Form("eta_%i_csv",j)]  -> Fill(jet->eta());
         h1[Form("phi_%i_csv",j)]  -> Fill(jet->phi());
         h1[Form("btag_%i_csv",j)] -> Fill(jet->btag());
      }
      if ( !isbbb ) h1["m12_csv"] -> Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
         
   }
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
// PRINT OUTS   
   
   // Cut flow
   // 0: triggered events
   // 1: 3+ idloose jets
   // 2: matched to online
   // 3: kinematics
   // 4: delta R
   // 5: delta eta
   // 6: btag (bbnb)
   
   double fracAbs[10];
   double fracRel[10];
   std::string cuts[10];
   cuts[0] = "HLT";
   cuts[1] = "Triple idloose-jet";
   cuts[2] = "Triple jet kinematics";
   cuts[3] = "Delta R(i;j)";
   cuts[4] = "Delta eta(j1;j2)";
   cuts[5] = "btagged (bbnb)";
   if ( isbbb ) cuts[5] = "btagged (bbb)";
   cuts[6] = "Matched to online j1;j2";
   
   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 7; ++i )
   {
      fracAbs[i] = double(nsel[i])/nsel[0];
      if ( i>0 )
         fracRel[i] = double(nsel[i])/nsel[i-1];
      else
         fracRel[i] = fracAbs[i];
      printf ("%-23s  %10d  %10.3f  %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
   }
   // CSV output
   printf ("%-23s , %10s , %10s , %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 7; ++i )
      printf ("%-23s , %10d , %10.3f , %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 

   // Trigger objects counts   
   std::cout << std::endl;
   printf ("%-40s  %10s \n", std::string("Trigger object").c_str(), std::string("# events").c_str() ); 
   for ( size_t io = 0; io < triggerObjects.size() ; ++io )
   {
      printf ("%-40s  %10d \n", triggerObjects[io].c_str(), nmatch[io] ); 
   }
   
   
   
   
      
//    
}

