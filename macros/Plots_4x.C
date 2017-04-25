void Plots_4x()
{
// // All had
//   TFile * f1 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw.root","old");
//   TFile * f2 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw.root","old");
// // All had MC
//   TFile * f3 = new TFile("qcd90x_2600_bx_HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw.root","old");
//   TFile * f4 = new TFile("qcd90x_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw.root","old");
   
// Semilep data
  TFile * f1 = new TFile("data_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw.root","old");
  TFile * f2 = new TFile("data_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw.root","old");
// Semilep MC
  TFile * f3 = new TFile("qcd90x_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw.root","old");
  TFile * f4 = new TFile("qcd90x_2600_bx_HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw.root","old");
   

   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   g1 -> SetName("rate1");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g2 -> SetName("rate2");
   TGraphAsymmErrors * g3 = (TGraphAsymmErrors*) f3 -> Get("rate_total");
   g3 -> SetName("rate3");
   TGraphAsymmErrors * g4 = (TGraphAsymmErrors*) f4 -> Get("rate_total");
   g4 -> SetName("rate4");
   
   int n = 2;
   double x[2] = {0,70};
   double y[2] = {0,0};
   TGraph * gdummy = new TGraph(n,x,y);
   
   g1 -> SetMarkerStyle(20);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   g3 -> SetMarkerStyle(25);
   g4 -> SetMarkerStyle(25);
   g4 -> SetMarkerColor(kRed);
   g4 -> SetLineColor(kRed);
   
   double minfit1 = 15;
   double maxfit1 = 40;
   double minfit2 = 30;
   double maxfit2 = 70;
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);
//   c1->SetLeftMargin(0.15);
   
   TF1 * pol1_1 = new TF1("pol1_1","pol1",minfit1,70);
   pol1_1 -> SetLineColor(kBlack);
   pol1_1 -> SetLineStyle(1);
   TFitResultPtr fit1 = g1 -> Fit(pol1_1,"S","",minfit1,maxfit1);
   
   TF1 * pol1_2 = new TF1("pol1_2","pol1",minfit1,70);
   pol1_2 -> SetLineColor(kRed);
   pol1_2 -> SetLineStyle(1);
   TFitResultPtr fit2 = g2 -> Fit(pol1_2,"S","",minfit1,maxfit1);
   
   
   TF1 * pol1_3 = new TF1("pol1_3","pol2",minfit2,70);
   pol1_3 -> SetLineColor(kBlack);
   pol1_3 -> SetLineStyle(3);
   TFitResultPtr fit3 = g3 -> Fit(pol1_3,"S","",minfit2,maxfit2);
   
   TF1 * pol1_4 = new TF1("pol1_4","pol1",minfit2,70);
   pol1_4 -> SetLineColor(kRed);
   pol1_4 -> SetLineStyle(3);
   TFitResultPtr fit4 = g4 -> Fit(pol1_4,"S","",minfit2,maxfit2);
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> Add(g1);
   mg -> Add(g2);
   mg -> Add(g3);
   mg -> Add(g4);
   mg -> Add(gdummy);
   
   
   mg -> Draw("ap");
   
   pol1_1 -> Draw("same");
   pol1_2 -> Draw("same");
   pol1_3 -> Draw("same");
   pol1_4 -> Draw("same");
   
   TLegend * leg = new TLegend(0.15,0.7,0.6,0.9); 
   leg->AddEntry(g1,"Data 2016 hltCSV > 0.84","p");
   leg->AddEntry(g2,"Data 2016 hltCSV > 0.92","p");
   leg->AddEntry(g3,"QCD 90X hltCSV > 0.84","p");
   leg->AddEntry(g4,"QCD 90X hltCSV > 0.92","p");
   
   leg -> Draw();
   mg -> GetXaxis() -> SetRangeUser(10,70);
   mg -> GetYaxis() -> SetTitle("rate (Hz)");
   mg -> GetXaxis() -> SetTitle("pile up");
//   mg -> GetYaxis() -> SetTitleOffset(2); 
   
   mg -> SetMinimum(0);
   mg -> SetMaximum(90);
   gPad->Modified();
   
   c1 -> SaveAs("data2016_qcd90x_rates_semilep.png");
   
}
