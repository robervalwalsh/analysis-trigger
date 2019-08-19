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
void reweightQCDtoDataCDE(OnlineBtagAnalyser &);
void reweightQCDtoDataF(OnlineBtagAnalyser &);

// =============================================================================================   
int main(int argc, char ** argv)
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   OnlineBtagAnalyser onlinebtag(argc,argv);
 
   // check if it is MC
   bool isMC = onlinebtag.config()->isMC();
   
   // eras and corresponding lumis
   auto eras = onlinebtag.config()->eras();
   auto eraslumi = onlinebtag.config()->erasLumi();
   
   // prescaling data
   int seedps = 4321;
   if ( onlinebtag.seed() <= 0 ) seedps += 1;
   else seedps += onlinebtag.seed();
   TRandom3 * rndps = nullptr;
   float ps = -1;
   if ( ! isMC && eras.size() > 0 )
   {
      if ( eras[0] == "2017C" and onlinebtag.config()->index() == 1 )
      {
         ps = 20;
         rndps = new TRandom3(seedps);
      }
   }
      
   int seed = 1234;
   if ( onlinebtag.seed() <= 0 ) seed += 1;
   else seed += onlinebtag.seed();
   
   
   // init for splitting MC into eras
   bool era_split = false; 
   float total_lumi = 0; 
   float fraction_accepted = 1; 
   TFile * fpu[2] = {nullptr,nullptr};
   TH1D  * hpu[2] = {nullptr,nullptr};
   TRandom3 * rnd = nullptr;
   if ( isMC && eras.size() == 2 && eraslumi.size() == 2 )
   {
      era_split = true;
      rnd = new TRandom3(seed);
      // first era is the nominal era, the otheris the rejected era
      total_lumi = eraslumi[0] + eraslumi[1];
      fraction_accepted = eraslumi[0]/total_lumi;
      //  fraction CDE / F
      total_lumi = 1.;
      if ( eras[0] == "2017CDE" ) fraction_accepted = 0.65;
      else fraction_accepted = 0.35;
      
      fpu[0] = new TFile(Form("/nfs/dust/cms/user/walsh/cms/analysis/cmssw/mssmhbb/2017/CMSSW_9_4_13/src/Analysis/Trigger/test/online-btag-2017_v2/qcd/MyDataPileupHistogram_Run%s.root",eras[0].c_str()));
      hpu[0] = (TH1D*) fpu[0] -> Get("pileup");
      hpu[0] -> SetName("pileup_0");
      hpu[0] -> Scale(1./hpu[0]->Integral());
      fpu[1] = new TFile(Form("/nfs/dust/cms/user/walsh/cms/analysis/cmssw/mssmhbb/2017/CMSSW_9_4_13/src/Analysis/Trigger/test/online-btag-2017_v2/qcd/MyDataPileupHistogram_Run%s.root",eras[1].c_str()));
      hpu[1] = (TH1D*) fpu[1] -> Get("pileup");
      hpu[1] -> SetName("pileup_1");
      hpu[1] -> Scale(1./hpu[1]->Integral());
   }
   
   onlinebtag.jetHistograms(2,"tag");
   onlinebtag.jetHistograms(2,"probe");
   onlinebtag.jetHistograms(2,"probe_match");
//   if ( isMC ) onlinebtag.pileupHistogram();
   
   onlinebtag.tagAndProbeTree();
   
   // Analysis of events
   std::cout << "The sample size is " << onlinebtag.analysis()->size() << " events" << std::endl;
   std::cout << "---------------------------" << std::endl;
   
   std::cout << "Workflow index = " << onlinebtag.config()->workflow() << std::endl;
   std::cout << "--------------------" << std::endl;
   
// 
   for ( int i = 0 ; i < onlinebtag.nEvents() ; ++i )
   {
      bool goodEvent = onlinebtag.event(i);

      if ( ! goodEvent ) continue;
      
   // splitting sample according to lumi and pileup profile
      if ( era_split )
      {
         auto hcut = onlinebtag.histogram("cutflow");
         auto cut = onlinebtag.cutflow();
         ++cut;
         if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
            hcut -> GetXaxis()-> SetBinLabel(cut+1,Form("Selecting according to era %s",eras[0].c_str()));
   
         float truepu = onlinebtag.trueInteractions();
         float pufrac[2];
         pufrac[0] = fraction_accepted      * hpu[0] -> GetBinContent(hpu[0] -> FindBin(truepu));
         pufrac[1] = (1.-fraction_accepted) * hpu[1] -> GetBinContent(hpu[1] -> FindBin(truepu));
         pufrac[0] /= (pufrac[0]+pufrac[1]);
         pufrac[1] /= (pufrac[0]+pufrac[1]);
         auto x = rnd->Rndm();
         
         if ( eras[0] == "2017CDE" ) 
         {
            if ( x >  pufrac[0] ) continue;
         }
         else
         {
            if ( x <  (1-pufrac[0]) ) continue;
         }
         
         hcut -> Fill(cut,onlinebtag.weight());
         onlinebtag.cutflow(cut);
//         onlinebtag.fillPileupHistogram();

      }

      
   // pileup weight
//      onlinebtag.actionApplyPileupWeight();
      
   // trigger selection
      if ( ! onlinebtag.selectionHLT()                 )   continue;
      if ( ! onlinebtag.selectionL1 ()                 )   continue;  // to be used in case of "OR" of seeds
   
      // prescale 2017C bin1
      if ( rndps )
      {
         auto hcut = onlinebtag.histogram("cutflow");
         auto cut = onlinebtag.cutflow();
         ++cut;
         if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
            hcut -> GetXaxis()-> SetBinLabel(cut+1,Form("Prescaling to %3.1f for era %s, bin %d",ps,eras[0].c_str(),onlinebtag.config()->index()));
         
         auto x = rndps->Rndm();
         if ( x > 1./ps ) continue;
         
         hcut -> Fill(cut,onlinebtag.weight());
         onlinebtag.cutflow(cut);
      }

         
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
   if ( eras[0] == "2017CDE" )  reweightQCDtoDataCDE(onlinebtag);
   if ( eras[0] == "2017F" )    reweightQCDtoDataF(onlinebtag);
      
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
   
   // this is not a good idea - use scale instead
//    if ( isMC && onlinebtag.config()->crossSection() > 0 )
//    {
//       float ngen = onlinebtag.histogram("cutflow") -> GetBinContent(4);
//       float xsection = onlinebtag.config()->crossSection();
//       float genlumi = ngen/xsection;
//       float scale = eraslumi[0]/genlumi;
//       onlinebtag.scale(scale);
//    }
   
} //end main



// ===========================
void reweightQCDtoDataF(OnlineBtagAnalyser & ana)
{
   if ( ! ana.config()->isMC() ) return;
   if ( ana.config()->index() < 1 ) return;
   
   auto hcut = ana.histogram("cutflow");
   auto cut = ana.cutflow();
   ++cut;
   if ( std::string(hcut -> GetXaxis()-> GetBinLabel(cut+1)) == "" )
      hcut -> GetXaxis()-> SetBinLabel(cut+1,"*** Re-weight QCD MC to data 2017F");
   
   // reweights for 2017 era F to match 2017 eras C+D+E
   TF1 wbin[6][6];
   // bin1 §
   wbin[0][0] = TF1("w_bin1_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[0][1] = TF1("w_bin1_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[0][2] = TF1("w_bin1_mb" ,"pol3(0)",0,1500);  // eta MB
   wbin[0][3] = TF1("w_bin1_pb" ,"pol3(0)",0,1500);  // eta PB
   wbin[0][4] = TF1("w_bin1_mbe","pol3(0)",0,1500);  // eta MBE
   wbin[0][5] = TF1("w_bin1_pbe","pol3(0)",0,1500);  // eta PBE
   // bin2  §
   wbin[1][0] = TF1("w_bin2_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[1][1] = TF1("w_bin2_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[1][2] = TF1("w_bin2_mb" ,"pol1(0)",0,1500);  // eta MB
   wbin[1][3] = TF1("w_bin2_pb" ,"pol1(0)",0,1500);  // eta PB
   wbin[1][4] = TF1("w_bin2_mbe","pol1(0)",0,1500);  // eta MBE
   wbin[1][5] = TF1("w_bin2_pbe","pol2(0)",0,1500);  // eta PBE
   // bin3  §
   wbin[2][0] = TF1("w_bin3_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[2][1] = TF1("w_bin3_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[2][2] = TF1("w_bin3_mb" ,"pol3(0)",0,1500);  // eta MB
   wbin[2][3] = TF1("w_bin3_pb" ,"pol3(0)",0,1500);  // eta PB
   wbin[2][4] = TF1("w_bin3_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[2][5] = TF1("w_bin3_pbe","pol2(0)",0,1500);  // eta PBE
   // bin4  §
   wbin[3][0] = TF1("w_bin4_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[3][1] = TF1("w_bin4_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[3][2] = TF1("w_bin4_mb" ,"pol3(0)",0,1500);  // eta MB
   wbin[3][3] = TF1("w_bin4_pb" ,"pol3(0)",0,1500);  // eta PB
   wbin[3][4] = TF1("w_bin4_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[3][5] = TF1("w_bin4_pbe","pol2(0)",0,1500);  // eta PBE
   // bin5  §
   wbin[4][0] = TF1("w_bin5_me" ,"pol1(0)",0,1500);  // eta ME
   wbin[4][1] = TF1("w_bin5_pe" ,"pol1(0)",0,1500);  // eta PE
   wbin[4][2] = TF1("w_bin5_mb" ,"pol1(0)",0,1500);  // eta MB
   wbin[4][3] = TF1("w_bin5_pb" ,"pol1(0)",0,1500);  // eta PB
   wbin[4][4] = TF1("w_bin5_mbe","pol1(0)",0,1500);  // eta MBE
   wbin[4][5] = TF1("w_bin5_pbe","pol2(0)",0,1500);  // eta PBE
   
   
   // parameters
   // bin 1
   // ME
   wbin[0][0].SetParameter(0,     -1.38873);
   wbin[0][0].SetParameter(1,    0.0550498);
   wbin[0][0].SetParameter(2, -0.000259024);
   // PE                      
   wbin[0][1].SetParameter(0,    -0.710186);
   wbin[0][1].SetParameter(1,    0.0357005);
   wbin[0][1].SetParameter(2, -0.000134505);
   // MB                      
   wbin[0][2].SetParameter(0,     -1.61038);
   wbin[0][2].SetParameter(1,    0.0802513);
   wbin[0][2].SetParameter(2, -0.000741898);
   wbin[0][2].SetParameter(3,  2.30218e-06);
   // PB                      
   wbin[0][3].SetParameter(0,     -1.65898);
   wbin[0][3].SetParameter(1,    0.0808882);
   wbin[0][3].SetParameter(2, -0.000744568);
   wbin[0][3].SetParameter(3,  2.33914e-06);
   // MBE                     
   wbin[0][4].SetParameter(0,     -1.68883);
   wbin[0][4].SetParameter(1,    0.0803238);
   wbin[0][4].SetParameter(2, -0.000709856);
   wbin[0][4].SetParameter(3,  2.09475e-06);
   // PBE                     
   wbin[0][5].SetParameter(0,     0.451228);
   wbin[0][5].SetParameter(1,  -0.00432294);
   wbin[0][5].SetParameter(2,  0.000330765);
   wbin[0][5].SetParameter(3, -1.91195e-06);
   
   // bin 2
   // ME
   wbin[1][0].SetParameter(0,     -5.77745);
   wbin[1][0].SetParameter(1,    0.0880766);
   wbin[1][0].SetParameter(2, -0.000279371);
   // PE                      
   wbin[1][1].SetParameter(0,     -6.69575);
   wbin[1][1].SetParameter(1,    0.0939525);
   wbin[1][1].SetParameter(2, -0.000277546);
   // MB                      
   wbin[1][2].SetParameter(0,     0.933951);
   wbin[1][2].SetParameter(1,  0.000433095);
   // PB                      
   wbin[1][3].SetParameter(0,     0.787343);
   wbin[1][3].SetParameter(1,   0.00138429);
   // MBE                     
   wbin[1][4].SetParameter(0,     0.551883);
   wbin[1][4].SetParameter(1,   0.00290096);
   // PBE                     
   wbin[1][5].SetParameter(0,     -3.57747);
   wbin[1][5].SetParameter(1,    0.0522595);
   wbin[1][5].SetParameter(2, -0.000141476);

   
   // bin 3
   // ME
   wbin[2][0].SetParameter(0,     -8.54811);    
   wbin[2][0].SetParameter(1,    0.0657555);    
   wbin[2][0].SetParameter(2, -0.000108932);   
   // PE
   wbin[2][1].SetParameter(0,     -8.56109);    
   wbin[2][1].SetParameter(1,    0.0675469);    
   wbin[2][1].SetParameter(2, -0.000115469);   
   // MB
   wbin[2][2].SetParameter(0,     -14.6869);    
   wbin[2][2].SetParameter(1,     0.144514);    
   wbin[2][2].SetParameter(2,  -0.00042885);    
   wbin[2][2].SetParameter(3,   4.1285e-07);    
   // PB
   wbin[2][3].SetParameter(0,     -19.4593);    
   wbin[2][3].SetParameter(1,     0.207611);    
   wbin[2][3].SetParameter(2, -0.000696133);    
   wbin[2][3].SetParameter(3,   7.7616e-07);    
   // MBE
   wbin[2][4].SetParameter(0,     -6.33536);    
   wbin[2][4].SetParameter(1,    0.0508149);    
   wbin[2][4].SetParameter(2, -8.48972e-05);    
   // PBE
   wbin[2][5].SetParameter(0,     -1.22657);    
   wbin[2][5].SetParameter(1,     0.014491);    
   wbin[2][5].SetParameter(2, -2.23136e-05);    
   
   // bin 4
   // ME
   wbin[3][0].SetParameter(0,     -17.8586);    
   wbin[3][0].SetParameter(1,    0.0821032);    
   wbin[3][0].SetParameter(2, -8.75879e-05);    
   // PE
   wbin[3][1].SetParameter(0,     -12.2956);    
   wbin[3][1].SetParameter(1,    0.0592316);    
   wbin[3][1].SetParameter(2, -6.51025e-05);   
   // MB
   wbin[3][2].SetParameter(0,     -63.9487);    
   wbin[3][2].SetParameter(1,     0.415683);    
   wbin[3][2].SetParameter(2, -0.000877524);    
   wbin[3][2].SetParameter(3,  6.11425e-07);    
   // PB
   wbin[3][3].SetParameter(0,     -45.4686);    
   wbin[3][3].SetParameter(1,     0.296697);    
   wbin[3][3].SetParameter(2, -0.000625097);    
   wbin[3][3].SetParameter(3,  4.34747e-07);    
   // MBE
   wbin[3][4].SetParameter(0,     -13.2135);    
   wbin[3][4].SetParameter(1,    0.0600719);    
   wbin[3][4].SetParameter(2, -6.19879e-05);    
   // PBE
   wbin[3][5].SetParameter(0,     -6.78826);    
   wbin[3][5].SetParameter(1,    0.0312724);    
   wbin[3][5].SetParameter(2, -3.01266e-05);    
    
   // bin 5
   // ME
   wbin[4][0].SetParameter(0,      2.87494);    
   wbin[4][0].SetParameter(1,  -0.00335089);    
   // PE
   wbin[4][1].SetParameter(0,     0.814192);    
   wbin[4][1].SetParameter(1,  5.22367e-05);    
   // MB
   wbin[4][2].SetParameter(0,      1.37201);    
   wbin[4][2].SetParameter(1, -0.000662398);    
   // PB
   wbin[4][3].SetParameter(0,      1.66209);    
   wbin[4][3].SetParameter(1,  -0.00118165);    
   // MBE
   wbin[4][4].SetParameter(0,      1.62254);    
   wbin[4][4].SetParameter(1,   -0.0010949);    
   // PBE
   wbin[4][5].SetParameter(0,     -1.41377);    
   wbin[4][5].SetParameter(1,   0.00865488);    
   wbin[4][5].SetParameter(2, -7.59207e-06);    
    
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


// ================================
void reweightQCDtoDataCDE(OnlineBtagAnalyser & ana)
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
   // bin1 §
   wbin[0][0] = TF1("w_bin1_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[0][1] = TF1("w_bin1_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[0][2] = TF1("w_bin1_mb" ,"pol2(0)",0,1500);  // eta MB
   wbin[0][3] = TF1("w_bin1_pb" ,"pol2(0)",0,1500);  // eta PB
   wbin[0][4] = TF1("w_bin1_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[0][5] = TF1("w_bin1_pbe","pol2(0)",0,1500);  // eta PBE
   // bin2  §
   wbin[1][0] = TF1("w_bin2_me" ,"pol2(0)",0,1500);  // eta ME
   wbin[1][1] = TF1("w_bin2_pe" ,"pol3(0)",0,1500);  // eta PE
   wbin[1][2] = TF1("w_bin2_mb" ,"pol2(0)",0,1500);  // eta MB
   wbin[1][3] = TF1("w_bin2_pb" ,"pol2(0)",0,1500);  // eta PB
   wbin[1][4] = TF1("w_bin2_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[1][5] = TF1("w_bin2_pbe","pol2(0)",0,1500);  // eta PBE
   // bin3  §
   wbin[2][0] = TF1("w_bin3_me" ,"pol3(0)",0,1500);  // eta ME
   wbin[2][1] = TF1("w_bin3_pe" ,"pol3(0)",0,1500);  // eta PE
   wbin[2][2] = TF1("w_bin3_mb" ,"pol3(0)",0,1500);  // eta MB
   wbin[2][3] = TF1("w_bin3_pb" ,"pol3(0)",0,1500);  // eta PB
   wbin[2][4] = TF1("w_bin3_mbe","pol2(0)",0,1500);  // eta MBE
   wbin[2][5] = TF1("w_bin3_pbe","expo+pol2(2)",0,1500);  // eta PBE
   // bin4  §
   wbin[3][0] = TF1("w_bin4_me" ,"pol1(0)",0,1500);  // eta ME
   wbin[3][1] = TF1("w_bin4_pe" ,"pol2(0)",0,1500);  // eta PE
   wbin[3][2] = TF1("w_bin4_mb" ,"pol2(0)",0,1500);  // eta MB
   wbin[3][3] = TF1("w_bin4_pb" ,"pol2(0)",0,1500);  // eta PB
   wbin[3][4] = TF1("w_bin4_mbe","pol1(0)",0,1500);  // eta MBE
   wbin[3][5] = TF1("w_bin4_pbe","pol3(0)",0,1500);  // eta PBE
   // bin5  §
   wbin[4][0] = TF1("w_bin5_me" ,"pol1(0)",0,1500);  // eta ME
   wbin[4][1] = TF1("w_bin5_pe" ,"pol1(0)",0,1500);  // eta PE
   wbin[4][2] = TF1("w_bin5_mb" ,"pol1(0)",0,1500);  // eta MB
   wbin[4][3] = TF1("w_bin5_pb" ,"pol1(0)",0,1500);  // eta PB
   wbin[4][4] = TF1("w_bin5_mbe","pol3(0)",0,1500);  // eta MBE
   wbin[4][5] = TF1("w_bin5_pbe","pol3(0)",0,1500);  // eta PBE
   
   
   // parameters
   // bin 1
   // ME
   wbin[0][0].SetParameter(0,-0.0652574);    
   wbin[0][0].SetParameter(1, 0.0229375);    
   wbin[0][0].SetParameter(2,-9.23584e-05);   //
   // PE
   wbin[0][1].SetParameter(0,-0.50297);    
   wbin[0][1].SetParameter(1, 0.0342216);    
   wbin[0][1].SetParameter(2,-0.000155978);   //
   // MB
   wbin[0][2].SetParameter(0, 0.257993);    
   wbin[0][2].SetParameter(1, 0.0165364);    
   wbin[0][2].SetParameter(2,-6.99772e-05);   //
   // PB
   wbin[0][3].SetParameter(0, 0.282757);    
   wbin[0][3].SetParameter(1, 0.0160886);    
   wbin[0][3].SetParameter(2,-6.94073e-05);  // 
   // MBE
   wbin[0][4].SetParameter(0, 0.0153423);    
   wbin[0][4].SetParameter(1, 0.0232781);    
   wbin[0][4].SetParameter(2,-0.000110542);   //
   // PBE
   wbin[0][5].SetParameter(0,-0.0344957);    
   wbin[0][5].SetParameter(1, 0.022747);    
   wbin[0][5].SetParameter(2,-9.22252e-05);   //
   
   // bin 2
   // ME
   wbin[1][0].SetParameter(0,-2.06559);    
   wbin[1][0].SetParameter(1, 0.0397254);    
   wbin[1][0].SetParameter(2,-0.000125515);  // 
   // PE
   wbin[1][1].SetParameter(0,-27.7487);    
   wbin[1][1].SetParameter(1, 0.565529);    
   wbin[1][1].SetParameter(2,-0.00364953);   
   wbin[1][1].SetParameter(3, 7.72979e-06);   //
   // MB
   wbin[1][2].SetParameter(0,-1.42803);    
   wbin[1][2].SetParameter(1, 0.0313777);    
   wbin[1][2].SetParameter(2,-9.909e-05);   //
   // PB
   wbin[1][3].SetParameter(0,-0.236463);    
   wbin[1][3].SetParameter(1, 0.0159934);    
   wbin[1][3].SetParameter(2,-5.05555e-05);   //
   // MBE
   wbin[1][4].SetParameter(0, 1.47188);    
   wbin[1][4].SetParameter(1,-0.00388074);    
   wbin[1][4].SetParameter(2, 4.9006e-06);   //
   // PBE
   wbin[1][5].SetParameter(0,-0.820059);    
   wbin[1][5].SetParameter(1, 0.0230815);    
   wbin[1][5].SetParameter(2,-7.16513e-05);  //
   
   // bin 3
   // ME
   wbin[2][0].SetParameter(0,     5.73088);    
   wbin[2][0].SetParameter(1,  -0.0642435);    
   wbin[2][0].SetParameter(2, 0.000271238);   
   wbin[2][0].SetParameter(3, -3.5903e-07);   
   // PE
   wbin[2][1].SetParameter(0,     -20.4788);    
   wbin[2][1].SetParameter(1,     0.213145);    
   wbin[2][1].SetParameter(2, -0.000686299);   
   wbin[2][1].SetParameter(3,  7.17865e-07);   
   // MB
   wbin[2][2].SetParameter(0,     -4.63744);    
   wbin[2][2].SetParameter(1,    0.0563069);    
   wbin[2][2].SetParameter(2, -0.000181273);    
   wbin[2][2].SetParameter(3,  1.87502e-07);    
   // PB
   wbin[2][3].SetParameter(0,     -9.06343);    
   wbin[2][3].SetParameter(1,     0.102401);    
   wbin[2][3].SetParameter(2, -0.000339927);    
   wbin[2][3].SetParameter(3,  3.68722e-07);    
   // MBE
   wbin[2][4].SetParameter(0,    0.580237);    
   wbin[2][4].SetParameter(1,  0.00263568);    
   wbin[2][4].SetParameter(2, -3.8462e-06);    
   // PBE
   wbin[2][5].SetParameter(0,    -0.834244);    
   wbin[2][5].SetParameter(1,   0.00802932);    
   wbin[2][5].SetParameter(2,     -7.39069);    
   wbin[2][5].SetParameter(3,    0.0648124);    
   wbin[2][5].SetParameter(4, -0.000175776);    
   
   // bin 4
   // ME
   wbin[3][0].SetParameter(0,    0.962966);    
   wbin[3][0].SetParameter(1, 1.66378e-05);    
   // PE
   wbin[3][1].SetParameter(0,     -7.81981);    
   wbin[3][1].SetParameter(1,    0.0396127);    
   wbin[3][1].SetParameter(2, -4.39654e-05);   
   // MB
   wbin[3][2].SetParameter(0,     -1.75287);    
   wbin[3][2].SetParameter(1,    0.0122508);    
   wbin[3][2].SetParameter(2, -1.34238e-05);    
   // PB
   wbin[3][3].SetParameter(0,     -2.78444);    
   wbin[3][3].SetParameter(1,    0.0172846);    
   wbin[3][3].SetParameter(2, -1.94909e-05);    
   // MBE
   wbin[3][4].SetParameter(0,     1.81482);    
   wbin[3][4].SetParameter(1, -0.00190951);    
   // PBE
   wbin[3][5].SetParameter(0,     -34.9409);    
   wbin[3][5].SetParameter(1,     0.233191);    
   wbin[3][5].SetParameter(2, -0.000499361);    
   wbin[3][5].SetParameter(3,  3.52979e-07);    
    
   // bin 5
   // ME
   wbin[4][0].SetParameter(0,     2.60688);    
   wbin[4][0].SetParameter(1, -0.00292853);    
   // PE
   wbin[4][1].SetParameter(0,      1.6879);    
   wbin[4][1].SetParameter(1, -0.00137564);    
   // MB
   wbin[4][2].SetParameter(0,      1.10502);    
   wbin[4][2].SetParameter(1, -0.000207952);    
   // PB
   wbin[4][3].SetParameter(0,      1.25404);    
   wbin[4][3].SetParameter(1, -0.000475007);    
   // MBE
   wbin[4][4].SetParameter(0,      22.8689);    
   wbin[4][4].SetParameter(1,    -0.107603);    
   wbin[4][4].SetParameter(2,  0.000175844);    
   wbin[4][4].SetParameter(3, -9.59547e-08);    
   // PBE
   wbin[4][5].SetParameter(0,      26.1834);    
   wbin[4][5].SetParameter(1,    -0.141815);    
   wbin[4][5].SetParameter(2,  0.000259597);    
   wbin[4][5].SetParameter(3, -1.54514e-07);    
    
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


