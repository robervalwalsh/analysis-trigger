#include "boost/program_options.hpp"
#include "boost/algorithm/string.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Core/interface/Analysis.h"
#include "Analysis/Core/bin/macro_config.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;

// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // read config parameters   
   if ( macro_config(argc, argv) != 0 ) return -1;
   
   // some specific requirements
   if ( njetsmin_ != 2 )
   {
      std::cout << "*** This is a tag and probe between two jets, nJetsMin must be equal to two. ***" << std::endl;
      return -1;
   }
   
   // Creat Analysis object
   Analysis analysis(inputlist_);
   
   // Jets
   analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsPuppi");
   // Trigger path info
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   // Trigger objects
   for ( auto & obj : triggerObjects_ )
   {
      analysis.addTree<TriggerObject> (obj,Form("MssmHbb/Events/slimmedPatTrigger/%s",obj.c_str()));
   }

   // JSON for data   
   if( !isMC_ ) analysis.processJsonFile(json_);
   
   // output file
   TFile hout(outputRoot_.c_str(),"recreate");
   
   // histograms
   std::map<std::string, TH1F*> h1;
   int nbins = 12+7+6+2;
   double ptbins[28] = {40,50,60,70,80,90,100,110,120,130,140,150,160,180,200,220,240,260,280,300,350,400,450,500,550,600,800,1000};
   
   h1["pt_tag"  ]      = new TH1F("pt_tag"      , "" , 100, 0   , 1000);
   h1["pt_tag_var"]    = new TH1F("pt_tag_var"  , "" , nbins, ptbins);
   h1["eta_tag" ]      = new TH1F("eta_tag"     , "" , 11 , -2.2, 2.2 );
   h1["phi_tag" ]      = new TH1F("phi_tag"     , "" , 16,  -3.2, 3.2 );
   h1["btag_tag"]      = new TH1F("btag_tag"    , "" , 200, 0   , 1   );
   h1["btaglog_tag"]   = new TH1F("btaglog_tag" , "" , 20 , 0   , 10  );
   
   h1["pt_probe_den"  ]    = new TH1F("pt_probe_den"      , "" , 100, 0   , 1000);
   h1["pt_probe_den_var"]  = new TH1F("pt_probe_den_var"  , "" , nbins, ptbins);
   h1["eta_probe_den" ]    = new TH1F("eta_probe_den"     , "" , 11 , -2.2, 2.2 );
   h1["phi_probe_den" ]    = new TH1F("phi_probe_den"     , "" , 16,  -3.2, 3.2 );
   h1["btag_probe_den"]    = new TH1F("btag_probe_den"    , "" , 200, 0   , 1   );
   h1["btaglog_probe_den"] = new TH1F("btaglog_probe_den" , "" , 20 , 0   , 10  );

   h1["pt_probe_num"  ]    = new TH1F("pt_probe_num"      , "" , 100, 0   , 1000);
   h1["pt_probe_num_var"]  = new TH1F("pt_probe_num_var"  , "" , nbins, ptbins);
   h1["eta_probe_num" ]    = new TH1F("eta_probe_num"     , "" , 11 , -2.2, 2.2 );
   h1["phi_probe_num" ]    = new TH1F("phi_probe_num"     , "" , 16,  -3.2, 3.2 );
   h1["btag_probe_num"]    = new TH1F("btag_probe_num"    , "" , 200, 0   , 1   );
   h1["btaglog_probe_num"] = new TH1F("btaglog_probe_num" , "" , 20 , 0   , 10  );
      
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
   
   int nsel[10] = { };
   int nmatch[10] = { };

   if ( nevtmax_ < 0 ) nevtmax_ = analysis.size();
   for ( int i = 0 ; i < nevtmax_ ; ++i )
   {
      bool goodEvent = true;
      
      if ( i > 0 && i%100000==0 ) std::cout << i << "  events processed! " << std::endl;
      
      analysis.event(i);
      if (! isMC_ )
      {
         if (!analysis.selectJson() ) continue; // To use only goodJSonFiles
      }
      
      int triggerFired = analysis.triggerResult(hltPath_);
      if ( !triggerFired ) continue;
      
//      std::cout << "oioi - trigger" << std::endl;
      
      ++nsel[0];
      
      // match offline to online
      analysis.match<Jet,TriggerObject>("Jets",triggerObjects_,0.5);
      
      // Jets - std::shared_ptr< Collection<Jet> >
      auto slimmedJets = analysis.collection<Jet>("Jets");
      std::vector<Jet *> selectedJets;
      for ( int j = 0 ; j < slimmedJets->size() ; ++j )
      {
         if ( jetsid_ == "TIGHT" ) // LOOSE is default if TIGHT not given
         {
            if ( slimmedJets->at(j).idTight() ) selectedJets.push_back(&slimmedJets->at(j));
         }
         else
         {
            if ( slimmedJets->at(j).idLoose() ) selectedJets.push_back(&slimmedJets->at(j));
         }
      }
      if ( (int)selectedJets.size() < njetsmin_ ) continue;
      
//      std::cout << "oioi - jet id" << std::endl;
      
      ++nsel[1];
      
      // kinematic and offline btagging selection - n leading jets
      for ( int j = 0; j < njetsmin_; ++j )
      {
         Jet * jet = selectedJets[j];
         if ( jet->pt() < jetsptmin_[j] || jet->pt() > jetsptmax_[j] || fabs(jet->eta()) > jetsetamax_[j] || jet->btag() < jetsbtagmin_[j] )
         {
            goodEvent = false;
            break;
         }
      }
      
      if ( ! goodEvent ) continue;
      
//      std::cout << "oioi - jet kinematics" << std::endl;
      
      ++nsel[2];
      
      // deltaR and deltaPhi
      for ( int j1 = 0; j1 < njetsmin_-1; ++j1 )
      {
         const Jet & jet1 = *selectedJets[j1];
         for ( int j2 = j1+1; j2 < njetsmin_; ++j2 )
         {
            const Jet & jet2 = *selectedJets[j2];
            if ( jet1.deltaR(jet2) < drmin_ ) goodEvent = false;
            if ( jet1.deltaPhi(jet2) < dphimin_ ) goodEvent = false;
         }
      }
      
      if ( ! goodEvent ) continue;
      
//      std::cout << "oioi - delta R and phi" << std::endl;
      
      ++nsel[3];
      
      
      // tag jets must be matched
      bool tagmatched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
      {       
         if ( ! selectedJets[1]->matched(triggerObjects_[io]) ) tagmatched[io] = false;
//         if ( matched[io] ) ++nmatch[io];
         goodEvent = ( goodEvent && tagmatched[io] );
      }
      
      if ( ! goodEvent ) continue;
      
//      std::cout << "oioi - tag match" << std::endl;
      
      ++nsel[4];
     
      // fill histograms for tag jet 
      h1["pt_tag"  ]    -> Fill(selectedJets[1]->pt()  );
      h1["pt_tag_var"]  -> Fill(selectedJets[1]->pt()  );
      h1["eta_tag" ]    -> Fill(selectedJets[1]->eta() );
      h1["phi_tag" ]    -> Fill(selectedJets[1]->phi() );
      h1["btag_tag"]    -> Fill(selectedJets[1]->btag());
      h1["btaglog_tag"] -> Fill(-log(1-selectedJets[1]->btag()));
      
      // fill histograms for probe jet not matched yet 
      h1["pt_probe_den"  ]    -> Fill(selectedJets[0]->pt()  );
      h1["pt_probe_den_var"]  -> Fill(selectedJets[0]->pt()  );
      h1["eta_probe_den" ]    -> Fill(selectedJets[0]->eta() );
      h1["phi_probe_den" ]    -> Fill(selectedJets[0]->phi() );
      h1["btag_probe_den"]    -> Fill(selectedJets[0]->btag());
      h1["btaglog_probe_den"] -> Fill(-log(1-selectedJets[0]->btag()));
      
      
      // is probe jet matched?
      bool probematched[10] = {true,true,true,true,true,true,true,true,true,true};
      for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
      {       
         if ( ! selectedJets[0]->matched(triggerObjects_[io]) ) probematched[io] = false;
         goodEvent = ( goodEvent && probematched[io] );
      }
      
      if ( ! goodEvent ) continue;
      
//      std::cout << "oioi - probe match" << std::endl;
      
      ++nsel[5];
      
      // fill histograms for probe jet 
      h1["pt_probe_num"  ]    -> Fill(selectedJets[0]->pt()  );
      h1["pt_probe_num_var"]  -> Fill(selectedJets[0]->pt()  );
      h1["eta_probe_num" ]    -> Fill(selectedJets[0]->eta() );
      h1["phi_probe_num" ]    -> Fill(selectedJets[0]->phi() );
      h1["btag_probe_num"]    -> Fill(selectedJets[0]->btag());
      h1["btaglog_probe_num"] -> Fill(-log(1-selectedJets[0]->btag()));
      
//            std::cout << "oioi - end" << std::endl;
//            std::cout << "oioi ============= " << std::endl;
      


   }
   
   for (auto & ih1 : h1)
   {
      ih1.second -> Write();
   }
   
   hout.Close();
   
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
//    if ( signalregion_ ) cuts[5] = "btagged (bbb)";
//    cuts[6] = "Matched to online j1;j2";
//    
   printf ("%-23s  %10s  %10s  %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
   for ( int i = 0; i < 6; ++i )
   {
      fracAbs[i] = double(nsel[i])/nsel[0];
      if ( i>0 )
         fracRel[i] = double(nsel[i])/nsel[i-1];
      else
         fracRel[i] = fracAbs[i];
      printf ("%-23s  %10d  %10.3f  %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
   }
//    // CSV output
//    printf ("%-23s , %10s , %10s , %10s \n", std::string("Cut flow").c_str(), std::string("# events").c_str(), std::string("absolute").c_str(), std::string("relative").c_str() ); 
//    for ( int i = 0; i < 7; ++i )
//       printf ("%-23s , %10d , %10.3f , %10.3f \n", cuts[i].c_str(), nsel[i], fracAbs[i], fracRel[i] ); 
// 
//    // Trigger objects counts   
//    std::cout << std::endl;
//    printf ("%-40s  %10s \n", std::string("Trigger object").c_str(), std::string("# events").c_str() ); 
// //   for ( size_t io = 0; io < triggerObjects_.size() ; ++io )
// //   {
// //      printf ("%-40s  %10d \n", triggerObjects_[io].c_str(), nmatch[io] ); 
// //   }
   
   
   
   
      
   
} //end main

