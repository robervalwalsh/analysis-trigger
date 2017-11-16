#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Tools/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


// =============================================================================================   
int main(int argc, char * argv[])
{
   bool isMC = false;
   bool isbbb = true;
   std::string inputList = "rootFilelist.txt";
   std::string outputRoot = "histograms_ntuples_purity_bbnb.root";
   std::string json = "json_2016H_r284044_ls1to30.txt";
   int bWP = 2;
   float btagcut[3] = {0.46,0.8434,0.935};
   // Cuts                                         // <<<<===== CMSDAS
   float ptmin[3]   = { 100.0, 100.0, 40.0 };
   float etamax[3]  = {   2.2,   2.2 , 2.2 };
   float btagmin[3] = { btagcut[bWP], btagcut[bWP], btagcut[bWP]};
   float nonbtag    = 0.46;
   float dRmin      = 1.;
   float detamax    = 1.55;
   float etacut_ = 2.3;
   float detacut_ = 1.6;
   
   
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   Analysis analysis(inputList,"MssmHbbTrigger/Events/EventInfo");
   
   analysis.addTree<Jet> ("Jets","MssmHbbTrigger/Events/slimmedJetsPuppi");
   analysis.addTree<Jet> ("hltJetsIdLoose","MssmHbbTrigger/Events/hltAK4PFJetsLooseIDCorrected");
   analysis.addTree<Jet> ("L1Jets","MssmHbbTrigger/Events/l1tJets");
   analysis.addTree<Jet> ("CaloJets","MssmHbbTrigger/Events/hltAK4CaloJetsCorrectedIDPassed");
   analysis.addTree<Jet> ("PFJets","MssmHbbTrigger/Events/hltAK4PFJetsCorrected");
   
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
   std::string hltPath = "HLT_DoubleJetsC100_DoubleBTagCSV_0p84_DoublePFJetsC100MaxDeta1p6_v";
   
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
   for ( int i = 0 ; i < analysis.size() ; ++i )
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
// L1 jets objects
      std::vector<Jet> l1jet100;
// L1 delta_eta jets objects
      std::vector<Jet> l1jet100deta;
// --------------------      
      // L1 jets selection
      for ( int j = 0; j < l1jets->size() ; ++j )
      {
         Jet l1j = l1jets->at(j);
         if ( l1j.pt() >= 100  && fabs(l1j.eta()) <= etacut_ ) l1jet100.push_back(l1j);
      }
      if ( l1jet100.size() < 2 ) continue;
      
      // L1 delta_eta jets selection
      if ( l1jet100.size() >= 2 )
      {
         for ( int j1 = 0 ; j1 < int(l1jet100.size())-1 ; ++j1 )
         {
            Jet jet1 = l1jet100.at(j1);
            for ( int j2 = j1+1 ; j2 < int(l1jet100.size()) ; ++j2 )
            {
               Jet jet2 = l1jet100.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= detacut_ )
               {
                  l1jet100deta.push_back(jet1);
                  l1jet100deta.push_back(jet2);
               }
            }
         }
      }
      if ( l1jet100deta.size() < 2 ) continue;
   
      
      // HLT
// --------------------      
      // Jets collections
      std::shared_ptr< Collection<Jet> >    calojets   = analysis.collection<Jet>("CaloJets");
      std::shared_ptr< Collection<Jet> >    pfjets     = analysis.collection<Jet>("PFJets");
      std::shared_ptr< Collection<JetTag> > jetstags   = analysis.collection<JetTag>("JetsTags");
// --------------------      
// HLT Jet objects
      // Calo
      std::vector<Jet> calojets100;
      // Calo for b-tagging
      std::vector<Jet> calojets80Selector;
      // PF
      std::vector<Jet> pfjets100;
      std::vector<Jet> pfjets100deta;
// HLT Jet btag      
      std::vector<JetTag>  bjets100wp092;
      std::vector<JetTag>  bjets100wp084;
      
      // -------------      
            
      // Jet selection
      for ( int j = 0 ; j < calojets->size() ; ++j )
      {
         Jet jet = calojets->at(j);
         if ( jet.pt() >= 100  && fabs(jet.eta()) <= etacut_ ) calojets100.push_back(jet);
         if ( jet.pt() >= 80  && fabs(jet.eta()) <= 2.4 ) calojets80Selector.push_back(jet);
      }
      if ( calojets100.size() < 2 ) continue;
      
      for ( int j = 0 ; j < pfjets->size() ; ++j )
      {
         Jet jet = pfjets->at(j);
         if ( jet.pt() >= 100  && fabs(jet.eta()) <= etacut_ ) pfjets100.push_back(jet);
      }
      if ( pfjets100.size() < 2 ) continue;
      
      
      // L1 delta_eta jets selection
      if ( pfjets100.size() >= 2 )
      {
         for ( int j1 = 0 ; j1 < int(pfjets100.size())-1 ; ++j1 )
         {
            Jet jet1 = pfjets100.at(j1);
            for ( int j2 = j1+1 ; j2 < int(pfjets100.size()) ; ++j2 )
            {
               Jet jet2 = pfjets100.at(j2);
               double deltaEta = fabs(jet1.eta() - jet2.eta());
               if ( deltaEta <= detacut_ )
               {
                  // semilep
                  pfjets100deta.push_back(jet1);
                  pfjets100deta.push_back(jet2);
               }
            }
         }
      }
      
      if ( pfjets100deta.size() < 2 ) continue;
      
      
      
      // online btagging
      std::vector<JetTag> bjets80wp092;
      for ( int bj = 0; bj < jetstags->size() ; ++bj )
      {
         JetTag bjet = jetstags->at(bj);
         if (bjet.pt() >= 80 && bj < 6 )
         {
            if ( bjet.btag() >= 0.92 ) bjets80wp092.push_back(bjet);
         }
      }
      if ( bjets80wp092.size() < 2 ) continue;
      
      
      // online jet id
      std::shared_ptr< Collection<Jet> > hltjetsidloose   = analysis.collection<Jet>("hltJetsIdLoose");
      std::vector<Jet> hltJetsIdLoose;
      for ( int j = 0; j < hltjetsidloose->size() ; ++j )
      {
         Jet jet = hltjetsidloose->at(j);
         if ( jet.pt() >= 100 && fabs(jet.eta()) <= 2.3 )
         {
            hltJetsIdLoose.push_back(jet);
         }
      }
//      if ( hltJetsIdLoose.size() < 2 ) continue;
      
      ++nsel[0];

      Vertex pv = offPVs->at(0);
      if ( !( ! pv.fake() && pv.ndof() > 4 && fabs(pv.z()) <= 24 && pv.rho() <= 2 ) ) continue;
      
      
      // match offline to online
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects,0.5);
      
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
      bool matched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( int j = 0; j < 2; ++j )
      {
         Jet * jet = selectedJets[j];
//         for ( auto & obj : triggerObjects )   matched = (matched && jet->matched(obj));
         for ( size_t io = 0; io < triggerObjects.size() ; ++io )
         {       
            if ( ! jet->matched(triggerObjects[io]) ) matched[io] = false;
         }
      }
      
      for ( size_t io = 0; io < triggerObjects.size() ; ++io )
      {
         if ( matched[io] ) ++nmatch[io];
         goodEvent = ( goodEvent && matched[io] );
      }
      
      if ( ! goodEvent ) continue;
      
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
   cuts[0] = "Triggered";
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

