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
#include "TF1.h" 

//#include "Analysis/Tools/interface/Analyser.h"
#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

using namespace std;
using namespace analysis;
//using namespace analysis::tools;
using namespace analysis::trigger;

using TH1s = std::map<std::string, TH1F*>;
using TH2s = std::map<std::string, TH2F*>;
     
void reweightFtoCDE(OnlineBtagAnalyser &);
void reweightQCDtoData(OnlineBtagAnalyser &);
void reweightQCDtoDataEta(OnlineBtagAnalyser &);

// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   OnlineBtagAnalyser onlinebtag(argc,argv);
   
   onlinebtag.jetHistograms(2,"tag");
   onlinebtag.jetHistograms(2,"probe");
   onlinebtag.jetHistograms(2,"probe_match");
   
   onlinebtag.tagAndProbeTree();
   
   // Analysis of events
   std::cout << "The sample size is " << onlinebtag.analysis()->size() << " events" << std::endl;
   std::cout << "---------------------------" << std::endl;
   
   std::cout << "Workflow index = " << onlinebtag.config()->workflow() << std::endl;
   std::cout << "--------------------" << std::endl;
   
   int seed = onlinebtag.seed();

// 
   for ( int i = 0 ; i < onlinebtag.nEvents() ; ++i )
   {
      bool goodEvent = onlinebtag.event(i);

      if ( ! goodEvent ) continue;

   // pileup weight
      onlinebtag.actionApplyPileupWeight();
      
   // trigger selection
      if ( ! onlinebtag.selectionHLT()                 )   continue;
      if ( ! onlinebtag.selectionL1 ()                 )   continue;  // to be used in case of "OR" of seeds
   
         
   // jet identification selection
      if ( ! onlinebtag.selectionJetId()               )   continue;
      if ( ! onlinebtag.selectionJetPileupId()         )   continue;
      if ( ! onlinebtag.selectionNJets()               )   continue;
      
      // jet energy resolution  
      onlinebtag.actionApplyJER();
      
   // delta Phi
      if ( ! onlinebtag.selectionJetDphi(1,2)          )   continue;
   // pT imbalance
      if ( ! onlinebtag.selectionJetPtImbalance(1,2)   )   continue;
      
   // tag jet selection
      if ( ! onlinebtag.selectionJet(2)                )   continue;
      if ( ! onlinebtag.onlineJetMatching(2)           )   continue;
      if ( ! onlinebtag.selectionBJet(2)               )   continue;
      // btag SF
      onlinebtag.actionApplyBtagSF(2);
      if ( ! onlinebtag.onlineBJetMatching(2)          )   continue;
   // additional tag selection
      if ( ! onlinebtag.selectionJetQGlikelihood(2)    )   continue;
      if ( ! onlinebtag.selectionBJetProbBB(2)         )   continue;
      if ( ! onlinebtag.selectionBJetProbC(2)          )   continue;
      
      onlinebtag.fillJetHistograms("tag");
      
   // probe jet selection
   // qg likelihood
      if ( ! onlinebtag.selectionJet(1)                )   continue;
      if ( ! onlinebtag.selectionJetQGlikelihood(1)    )   continue;
      if ( ! onlinebtag.onlineJetMatching(1)           )   continue;
      if ( ! onlinebtag.selectionBJet(1)               )   continue;
      // btag SF
      onlinebtag.actionApplyBtagSF(1);
      
   // reweight 2017F - too weird!?
//      reweightFtoCDE(onlinebtag);
   // reweight QCD to data CDE
      reweightQCDtoDataEta(onlinebtag);
      
   // If there is any AI selection, otherwise always true
      if ( ! onlinebtag.selectionAI()                  )   continue;
   // fill tree for AI
      onlinebtag.fillTagAndProbeTree();

   // fill tag histograms after selection
      onlinebtag.fillJetHistograms("probe");
      
   // PROBE jet match to online
      if ( ! onlinebtag.onlineBJetMatching(1)          )   continue;
      
   // fill probe histograms after selection
      onlinebtag.fillJetHistograms("probe_match");
   }
   
} //end main

void reweightQCDtoDataEta(OnlineBtagAnalyser & ana)
{
   if ( ! ana.config()->isMC() ) return;
   if ( ana.config()->index() < 1 ) return;
   
   auto hcut = ana.histogram("cutflow");
   auto cut = ana.cutflow();
   ++cut;
   if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
      hcut -> GetXaxis()-> SetBinLabel(cut+1,"*** Re-weight QCD MC to data 2017CDE");
   
   // reweights for 2017 era F to match 2017 eras C+D+E
   TF1 wbin[6][6];
   // bin1
   wbin[0][0] = TF1("w_bin1_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[0][1] = TF1("w_bin1_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[0][2] = TF1("w_bin1_mb" ,"pol2(0)",0,1500);  // eta MB
   wbin[0][3] = TF1("w_bin1_pb" ,"pol2(0)",0,1500);  // eta PB
   wbin[0][4] = TF1("w_bin1_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[0][5] = TF1("w_bin1_pbe","pol2(0)",0,1500);  // eta PBE
   // bin2
   wbin[1][0] = TF1("w_bin2_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[1][1] = TF1("w_bin2_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[1][2] = TF1("w_bin2_mb" ,"[0]",0,1500);  // eta MB
   wbin[1][3] = TF1("w_bin2_pb" ,"[0]",0,1500);  // eta PB
   wbin[1][4] = TF1("w_bin2_mbe","[0]",0,1500);  // eta MBE
   wbin[1][5] = TF1("w_bin2_pbe","pol2(0)",0,1500);  // eta PBE
   // bin3
   wbin[2][0] = TF1("w_bin3_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[2][1] = TF1("w_bin3_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[2][2] = TF1("w_bin3_mb" ,"[0]",0,1500);  // eta MB
   wbin[2][3] = TF1("w_bin3_pb" ,"[0]",0,1500);  // eta PB
   wbin[2][4] = TF1("w_bin3_mbe","[0]",0,1500);  // eta MBE
   wbin[2][5] = TF1("w_bin3_pbe","[0]",0,1500);  // eta PBE
   // bin4
   wbin[3][0] = TF1("w_bin4_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[3][1] = TF1("w_bin4_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[3][2] = TF1("w_bin4_mb" ,"[0]",0,1500);  // eta MB
   wbin[3][3] = TF1("w_bin4_pb" ,"[0]",0,1500);  // eta PB
   wbin[3][4] = TF1("w_bin4_mbe","pol1(0)",0,1500);  // eta MBE
   wbin[3][5] = TF1("w_bin4_pbe","[0]",0,1500);  // eta PBE
   // bin5
   wbin[4][0] = TF1("w_bin5_me" ,"pol1(0)",0,1500);  // eta ME
   wbin[4][1] = TF1("w_bin5_pe" ,"pol1(0)",0,1500);  // eta PE
   wbin[4][2] = TF1("w_bin5_mb" ,"pol1(0)",0,1500);  // eta MB
   wbin[4][3] = TF1("w_bin5_pb" ,"pol1(0)",0,1500);  // eta PB
   wbin[4][4] = TF1("w_bin5_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[4][5] = TF1("w_bin5_pbe","pol2(0)",0,1500);  // eta PBE
   // bin6
   wbin[5][0] = TF1("w_bin6_me" ,"pol1(0)",0,1500);  // eta ME
   wbin[5][1] = TF1("w_bin6_pe" ,"pol1(0)",0,1500);  // eta PE
   wbin[5][2] = TF1("w_bin6_mb" ,"pol1(0)",0,1500);  // eta MB
   wbin[5][3] = TF1("w_bin6_pb" ,"pol1(0)",0,1500);  // eta PB
   wbin[5][4] = TF1("w_bin6_mbe","pol1(0)",0,1500);  // eta MBE
   wbin[5][5] = TF1("w_bin6_pbe","pol1(0)",0,1500);  // eta PBE
   
   
   // parameters
   // bin 1
   // ME
   wbin[0][0].SetParameter(0,-2.37E-1);    
   wbin[0][0].SetParameter(1, 2.95E-2);    
   wbin[0][0].SetParameter(2,-1.49E-4);   
   // PE
   wbin[0][1].SetParameter(0,-1.90E-1);    
   wbin[0][1].SetParameter(1, 2.71E-2);    
   wbin[0][1].SetParameter(2,-1.25E-4);   
   // MB
   wbin[0][2].SetParameter(0, 4.15E-1);    
   wbin[0][2].SetParameter(1, 1.31E-2);    
   wbin[0][2].SetParameter(2,-5.71E-5);   
   // PB
   wbin[0][3].SetParameter(0, 3.29E-1);    
   wbin[0][3].SetParameter(1, 1.52E-2);    
   wbin[0][3].SetParameter(2,-6.66E-5);   
   // MBE
   wbin[0][4].SetParameter(0, 2.02E-1);    
   wbin[0][4].SetParameter(1, 1.80E-2);    
   wbin[0][4].SetParameter(2,-7.75E-5);   
   // PBE
   wbin[0][5].SetParameter(0, 2.60E-1);    
   wbin[0][5].SetParameter(1, 1.69E-2);    
   wbin[0][5].SetParameter(2,-7.61E-5);   
   
   // bin 2
   // ME
   wbin[1][0].SetParameter(0,-1.83E0);    
   wbin[1][0].SetParameter(1, 3.75E-2);    
   wbin[1][0].SetParameter(2,-1.21E-4);   
   // PE
   wbin[1][1].SetParameter(0,-2.00E0);    
   wbin[1][1].SetParameter(1, 3.88E-2);    
   wbin[1][1].SetParameter(2,-1.23E-4);   
   // MB
   wbin[1][2].SetParameter(0, 1.);    
   // PB
   wbin[1][3].SetParameter(0, 1.);    
   // MBE
   wbin[1][4].SetParameter(0, 1.);    
   // PBE
   wbin[1][5].SetParameter(0,-15.7E-1);    
   wbin[1][5].SetParameter(1, 3.28E-2);    
   wbin[1][5].SetParameter(2,-1.02E-4);  
   
   // bin 3
   // ME
   wbin[2][0].SetParameter(0,-2.84E0);    
   wbin[2][0].SetParameter(1, 2.68E-2);    
   wbin[2][0].SetParameter(2,-4.52E-5);   
   // PE
   wbin[2][1].SetParameter(0,-2.42E0);    
   wbin[2][1].SetParameter(1, 2.51E-2);    
   wbin[2][1].SetParameter(2,-4.48E-5);   
   // MB
   wbin[2][2].SetParameter(0, 1.);    
   // PB
   wbin[2][3].SetParameter(0, 1.);    
   // MBE
   wbin[2][4].SetParameter(0, 1.);    
   // PBE
   wbin[2][5].SetParameter(0, 1.);    
   
   // bin 4
   // ME
   wbin[3][0].SetParameter(0, 3.98E0);    
   wbin[3][0].SetParameter(1,-1.50E-2);    
   wbin[3][0].SetParameter(2, 1.85E-5);   
   // PE
   wbin[3][1].SetParameter(0,-3.84E0);    
   wbin[3][1].SetParameter(1, 2.13E-2);    
   wbin[3][1].SetParameter(2,-2.30E-5);   
   // MB
   wbin[3][2].SetParameter(0, 1.);    
   // PB
   wbin[3][3].SetParameter(0, 1.);    
   // MBE
   wbin[3][4].SetParameter(0, 1.51E0);    
   wbin[3][4].SetParameter(1,-1.19E-3);    
   // PBE
   wbin[3][5].SetParameter(0, 1.);    
    
   // bin 5
   // ME
   wbin[4][0].SetParameter(0, 2.66E0);    
   wbin[4][0].SetParameter(1,-2.98E-3);    
   // PE
   wbin[4][1].SetParameter(0, 1.74E0);    
   wbin[4][1].SetParameter(1,-1.37E-3);    
   // MB
   wbin[4][2].SetParameter(0, 1.29E0);    
   wbin[4][2].SetParameter(1,-5.28E-4);    
   // PB
   wbin[4][3].SetParameter(0, 1.45E0);    
   wbin[4][3].SetParameter(1,-8.20E-4);    
   // MBE
   wbin[4][4].SetParameter(0, 4.75E0);    
   wbin[4][4].SetParameter(1,-1.17E-2);    
   wbin[4][4].SetParameter(2, 8.77E-6);    
   // PBE
   wbin[4][5].SetParameter(0,-4.09E0);    
   wbin[4][5].SetParameter(1, 1.72E-2);    
   wbin[4][5].SetParameter(2,-1.42E-5);    
    
   // bin 6
   // ME
   wbin[5][0].SetParameter(0, 2.55E0);    
   wbin[5][0].SetParameter(1,-2.17E-3);    
   // PE
   wbin[5][1].SetParameter(0, 2.93E0);    
   wbin[5][1].SetParameter(1,-2.75E-3);    
   // MB
   wbin[5][2].SetParameter(0, 1.59E0);    
   wbin[5][2].SetParameter(1,-8.41E-4);    
   // PB
   wbin[5][3].SetParameter(0, 1.59E0);    
   wbin[5][3].SetParameter(1,-8.01E-4);    
   // MBE
   wbin[5][4].SetParameter(0, 1.26E0);    
   wbin[5][4].SetParameter(1,-3.72E-4);    
   // PBE
   wbin[5][5].SetParameter(0, 5.32E-1);    
   wbin[5][5].SetParameter(1, 7.13E-4);    
    
   // reweight QCD
   auto jets = ana.selectedJets();
   float pt  = jets[1]->pt();  // tag jet index = 1; probe jet index = 0
   float eta = jets[1]->eta();  // tag jet index = 1; probe jet index = 0
   int eta_indx = -1;
   if ( eta < -1.4 )               eta_indx = 0;
   if ( eta >  1.4 )               eta_indx = 1;
   if ( eta > -1.0 && eta <  0.0 ) eta_indx = 2; 
   if ( eta <  1.0 && eta >  0.0 ) eta_indx = 3; 
   if ( eta > -1.4 && eta < -1.0 ) eta_indx = 4;
   if ( eta <  1.4 && eta >  1.0 ) eta_indx = 5;

   int windex = ana.config()->index()-1;
   float myweight = 1;
   myweight = ana.weight()*(wbin[windex][eta_indx].Eval(pt));
   
   ana.weight(myweight);
   
   hcut -> Fill(cut,myweight);
   ana.cutflow(cut);
}



void reweightQCDtoData(OnlineBtagAnalyser & ana)
{
   if ( ! ana.config()->isMC() ) return;
   if ( ana.config()->index() < 1 ) return;
   
   auto hcut = ana.histogram("cutflow");
   auto cut = ana.cutflow();
   ++cut;
   if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
      hcut -> GetXaxis()-> SetBinLabel(cut+1,"*** Re-weight QCD MC to data 2017CDE");
   
   // reweights for 2017 era F to match 2017 eras C+D+E
   TF1 wbin[6];
   wbin[0] = TF1("w_bin1","pol2(0)",0,1500);
   wbin[1] = TF1("w_bin2","[0]",0,1500);
   wbin[2] = TF1("w_bin3","pol2(0)",0,1500);
   wbin[3] = TF1("w_bin4","[0]",0,1500);
   wbin[4] = TF1("w_bin5","[0]",0,1500);
   wbin[5] = TF1("w_bin6","expo(0)+pol0(2)",0,1500);
   
   // bin 1
   wbin[0].SetParameter(0,1.90E-1);    // +- 0.61E-1
   wbin[0].SetParameter(1,1.87E-2);    // +- 0.14E-2
   wbin[0].SetParameter(2,-8.61E-5);   // +- 0.80E-5
   // bin 2
   wbin[1].SetParameter(0,1.);
   // bin 3
   wbin[2].SetParameter(0,-7.48E-1);   // +- 3.23E-1
   wbin[2].SetParameter(1,1.26E-2);    // +- 0.23E-2
   wbin[2].SetParameter(2,-2.21E-5);   // +- 0.40E-5
   // bin 4
   wbin[3].SetParameter(0,1.);
   // bin 5
   wbin[4].SetParameter(0,1);
   // bin 6
   wbin[5].SetParameter(0,6.2);        // +- 10.3
   wbin[5].SetParameter(1,-1.44E-2);   // +- 2.09E-2
   wbin[5].SetParameter(2,9.19E-1);    // +- 0.98E-1
   
   // reweight QCD
   auto jets = ana.selectedJets();
   float jetpt = jets[1]->pt();  // tag jet index = 1; probe jet index = 0
   int windex = ana.config()->index()-1;
   float myweight = ana.weight()*(wbin[windex].Eval(jetpt));
   ana.weight(myweight);
   
   hcut -> Fill(cut,myweight);
   ana.cutflow(cut);
}

void reweightFtoCDE(OnlineBtagAnalyser & ana)
{
   if ( ana.config()->isMC() ) return;
   if ( ! ( ana.analysis()->run() >= 305040 && ana.analysis()->run() <= 306462 ) ) return;
   
   auto hcut = ana.histogram("cutflow");
   auto cut = ana.cutflow();
   ++cut;
   if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
      hcut -> GetXaxis()-> SetBinLabel(cut+1,"*** Re-weight data 2017F to data 2017CDE");
   
   // reweights for 2017 era F to match 2017 eras C+D+E
   TF1 wbin[6];
   wbin[0] = TF1("w_bin1","pol3(0)",0,1500);
   wbin[1] = TF1("w_bin2","pol1(0)",0,1500);
   wbin[2] = TF1("w_bin3","expo(0)+pol0(2)",0,1500);
   wbin[3] = TF1("w_bin4","expo(0)+pol0(2)",0,1500);
   wbin[4] = TF1("w_bin5","[0]",0,1500);
   wbin[5] = TF1("w_bin6","[0]",0,1500);
   
   // bin 1
   wbin[0].SetParameter(0,3.31);           // +- 0.53
   wbin[0].SetParameter(1,-7.23E-2);       // +- 1.90E-2
   wbin[0].SetParameter(2,6.98E-4);        // +- 2.16E-4
   wbin[0].SetParameter(3,-2.24E-6);       // +- 0.78E-6
   // bin 2                                
   wbin[1].SetParameter(0,1.405);          // +- 0.083
   wbin[1].SetParameter(1,-2.68E-3);       // +- 0.53E-3
   // bin 3                                
   wbin[2].SetParameter(0,8.9);            // +- 2.1
   wbin[2].SetParameter(1,-4.44E-2);       // +- 0.92E-2
   wbin[2].SetParameter(2,8.44E-2);        // +- 0.20E-2
   // bin 4                                
   wbin[3].SetParameter(0,21.7);           // +- 5.1
   wbin[3].SetParameter(1,-6.0E-2);        // +- 1.4E-2
   wbin[3].SetParameter(2,9.01E-1);        // +- 0.20E-1
   // bin 5                                 
   wbin[4].SetParameter(0,1);               
   // bin 6                                 
   wbin[5].SetParameter(0,1);
   
   // reweight 2017F
   auto jets = ana.selectedJets();
   float jetpt = jets[1]->pt();  // tag jet index = 1; probe jet index = 0
   int windex = ana.config()->workflow()-1;
   float myweight = ana.weight()*(wbin[windex].Eval(jetpt));
   ana.weight(myweight);
   
   hcut -> Fill(cut,myweight);
   ana.cutflow(cut);
   
}


