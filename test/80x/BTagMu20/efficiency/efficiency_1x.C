void efficiency_1x()
{
   TFile*  f = new TFile("mssmhbb_triggers_data_test.root","old");
   TH1F * h1 = (TH1F*) f -> Get("muon12");
   TH1F * h2 = (TH1F*) f -> Get("muon12Trg");
   
   h1 -> Rebin(2);
   h2 -> Rebin(2);
   
   TGraphAsymmErrors * g1 = new TGraphAsymmErrors(h2,h1,"cl=0.683 b(1,1) mode");
   
   g1 -> SetMarkerStyle(20);
   g1 -> SetMarkerColor(kBlack);
   g1 -> SetLineColor(kBlack);
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);

   TMultiGraph * mg = new TMultiGraph();
   mg -> Add(g1);
   
   mg -> Draw("AP");
   
   TLegend * leg = new TLegend(0.45,0.15,0.85,0.3);
   leg->SetHeader("online csv > 0.92");
   leg->AddEntry(g1,"jet pT > 40 GeV","pl");
   
//   leg -> Draw();
   mg -> GetYaxis() -> SetTitle("efficiency");
   mg -> GetXaxis() -> SetTitle("muon  pt(GeV)");
   mg -> GetXaxis() -> SetRangeUser(0,100);
   mg -> GetXaxis() -> SetTitleOffset(1.2); 
   mg -> GetYaxis() -> SetTitleOffset(1.3); 
   
   gPad->Modified();
   
   c1 -> SaveAs("muon12_eff_data80x_btagmu.png");
   
}
