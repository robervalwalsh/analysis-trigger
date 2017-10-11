#include "tdrstyle.C"

int EffsPlot()
{
   setTDRStyle();

   int min = 5;
   int max = 6;
   
   TFile * f[10];
   TGraphAsymmErrors * g_eff[10];
   
   f[0] = TFile::Open("histograms_2017C-v1_tep_40to1000_all.root");
   f[1] = TFile::Open("histograms_2017C-v2_tep_40to1000_all.root");
   f[2] = TFile::Open("histograms_2017C-v3_tep_40to1000_all.root");
   f[3] = TFile::Open("histograms_2017D-v1_tep_40to1000_all.root");
   
   f[5] = TFile::Open("histograms_2017C_tep_40to1000_all.root");
   f[6] = TFile::Open("histograms_2017D-v1_tep_40to1000_all.root");
   
   
   f[9] = TFile::Open("histograms_2017CD_tep_40to100_all.root");
   
   TCanvas * c1 = new TCanvas("c1","");
   
   bool first = true;
   int mycolor = 1;
   
   for ( int i = min; i < max; ++i )
   {
      TH1F * h_num = (TH1F*) f[i]->Get("pt_probe_num_var");
      TH1F * h_den = (TH1F*) f[i]->Get("pt_probe_den_var");
   
      g_eff[i] = new TGraphAsymmErrors(h_num,h_den,"cl=0.683 b(1,1) mode");
   
      if ( first )
      {
         g_eff[i] -> Draw("AP");
         g_eff[i] -> GetXaxis() -> SetTitle("probe jet p_{T} (GeV)");
//         g_eff[i] -> GetXaxis() -> SetTitle("probe jet eta");
//         g_eff[i] -> GetXaxis() -> SetTitle("probe jet -ln(1-CSVv2)");
         g_eff[i] -> GetYaxis() -> SetTitle("relative online efficiency");
         g_eff[i] -> GetYaxis() -> SetRangeUser(0,1);
         first = false;
      }
      else
      {
         g_eff[i] -> Draw("P");
      }
      if ( (max-min)==1 )
         mycolor = 1;
      
      if ( mycolor == 3 ) mycolor = 4;

      g_eff[i] -> SetMarkerStyle(20);
      g_eff[i] -> SetLineColor(mycolor);
      g_eff[i] -> SetMarkerColor(mycolor);
      
      mycolor += 1;
    
   }
   
   c1 -> SaveAs("eff_pt_2017C_40to100.png");
   
   return 0;
}
