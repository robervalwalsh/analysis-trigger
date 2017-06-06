void Plots_1xData()
{
   // All had
   TFile * f1 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw.root","old");
   TFile * f2 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw.root","old");
   
   // Semilep
//   TFile * f1 = new TFile("data_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw.root","old");
//   TFile * f2 = new TFile("data_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw.root","old");
   

   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   g1 -> SetName("rate1");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g2 -> SetName("rate2");
   
   int n = 2;
   double x[2] = {0,70};
   double y[2] = {0,0};
   TGraph * gdummy = new TGraph(n,x,y);
   
   g1 -> SetMarkerStyle(20);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);
//   c1->SetLeftMargin(0.15);
   
   TF1 * pol1_1 = new TF1("pol1_1","pol1",15,70);
   pol1_1 -> SetLineColor(kBlack);
   TFitResultPtr fit1 = g1 -> Fit(pol1_1,"S","",15,40);
   
   TF1 * pol1_2 = new TF1("pol1_2","pol1",15,70);
   pol1_2 -> SetLineColor(kRed);
   TFitResultPtr fit2 = g2 -> Fit(pol1_2,"S","",15,40);
   
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> Add(g1);
   mg -> Add(g2);
   mg -> Add(gdummy);
   
   
   mg -> Draw("ap");
   
   pol1_1 -> Draw("same");
   pol1_2 -> Draw("same");
   
   mg -> GetXaxis() -> SetRangeUser(10,70);
   mg -> GetYaxis() -> SetTitle("rate (Hz)");
   mg -> GetXaxis() -> SetTitle("pile up");
//   mg -> GetYaxis() -> SetTitleOffset(2); 
   
   mg -> SetMinimum(0);
   mg -> SetMaximum(40);
   gPad->Modified();
   
   c1 -> SaveAs("data2016_rates_allhad.png");
   
}
