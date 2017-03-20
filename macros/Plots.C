void Plots()
{
   TFile * f1 = new TFile("L1_DoubleJetC100.root","old");
   TFile * f2 = new TFile("L1_DoubleJetC100Eta2p3_dEtaMax1p6.root","old");
   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   TF1 * pol2_1 = new TF1("pol2_1","pol2",15,70);
   pol2_1 -> SetLineColor(kBlack);
   pol2_1 -> SetLineStyle(2);
   pol2_1 -> SetLineWidth(3);
   TFitResultPtr fit1 =	g1 -> Fit(pol2_1,"S");
   TF1 * pol2_2 = new TF1("pol2_2","pol2",15,70);
   pol2_2-> SetLineColor(kRed);
   pol2_2 -> SetLineStyle(2);
   pol2_2 -> SetLineWidth(3);
   TFitResultPtr fit2 =	g2 -> Fit(pol2_2,"S");
   
 
   TFile * fd1 = new TFile("data_L1_DoubleJetC100.root","old");
   TFile * fd2 = new TFile("data_L1_DoubleJetC100Eta2p3_dEtaMax1p6.root","old");
   TGraphAsymmErrors * gd1 = (TGraphAsymmErrors*) fd1 -> Get("rate_total");
   TGraphAsymmErrors * gd2 = (TGraphAsymmErrors*) fd2 -> Get("rate_total");
   
   gd2 -> SetMarkerColor(kRed);
   
 // WBM L1 functions   
   TF1 * func1 = new TF1("f1","217 + 98.7*x + 0.381*x*x",0,70);
   func1 -> SetLineColor(kBlack);
   func1 -> SetLineWidth(2);
  
   TCanvas * c1 = new TCanvas("c1","", 700,600);
   c1->SetLeftMargin(0.15);
   
   TH1F * h_c1 = c1->DrawFrame(15,0,70,20000);
   h_c1 -> GetYaxis() -> SetTitleOffset(2);
   h_c1 -> GetYaxis() -> SetTitle("rate (Hz)");
   h_c1 -> GetXaxis() -> SetTitle("pile up");
   
//   h_c1 -> Draw();
   gd1 -> Draw("p");
   gd2 -> Draw("p");
   
   pol2_1 -> Draw("same");
   pol2_2 -> Draw("same");
   
   func1 -> Draw("same");
   
   TLegend * leg = new TLegend(0.15,0.7,0.9,0.9);   
   leg->AddEntry(gd1,"L1_DoubleJetC100 - Parked ZeroBias","p");
   leg->AddEntry(pol2_1,"L1_DoubleJetC100 - QCD MC (fit pol2)","l");
   leg->AddEntry(func1,"L1_DoubleJetC100 - WBM data (fit pol2)","l");
   leg->AddEntry(gd2,"L1_DoubleJet100_Eta2p3_Deta1p6 - Parked ZeroBias","p");
   leg->AddEntry(pol2_2,"L1_DoubleJet100_Eta2p3_Deta1p6 - QCD MC (fit pol2)","l");
   
   
   leg -> Draw();
//   mg -> Draw("a");
//   mg -> GetXaxis() -> SetRangeUser(20,70);
//    mg -> SetMinimum(0);
//    mg -> SetMaximum(10000);
//    
    gPad->Modified();
 
}
