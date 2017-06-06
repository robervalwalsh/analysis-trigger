void efficiency()
{
   TFile*  f = new TFile("mssmhbb_triggers_TTbar.root","old");
   TH1F * h1 = (TH1F*) f -> Get("jet40mubtagTrg");
   TH1F * h2 = (TH1F*) f -> Get("jet40btagTrg");
   TH1F * h3 = (TH1F*) f -> Get("jet40mubtag");
   TH1F * h4 = (TH1F*) f -> Get("jet40btag");
   
   h1 -> Rebin(2);
   h2 -> Rebin(2);
   h3 -> Rebin(2);
   h4 -> Rebin(2);
   
   TGraphAsymmErrors * g1 = new TGraphAsymmErrors(h1,h3,"cl=0.683 b(1,1) mode");
   TGraphAsymmErrors * g2 = new TGraphAsymmErrors(h2,h4,"cl=0.683 b(1,1) mode");
   
   g1 -> SetMarkerStyle(20);
   g1 -> SetMarkerColor(kBlack);
   g1 -> SetLineColor(kBlack);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);

   TMultiGraph * mg = new TMultiGraph();
   mg -> Add(g1);
   mg -> Add(g2);
   
   mg -> Draw("AP");
   
   TLegend * leg = new TLegend(0.45,0.15,0.85,0.3);
   leg->SetHeader("online csv > 0.92");
   leg->AddEntry(g1,"semileptonic jet","pl");
   leg->AddEntry(g2,"hadronic jet","pl");
   
   leg -> Draw();
   mg -> GetYaxis() -> SetTitle("efficiency");
   mg -> GetXaxis() -> SetTitle("-ln(1-offline_CSVv2)");
   mg -> GetXaxis() -> SetTitleOffset(1.2); 
   mg -> GetYaxis() -> SetTitleOffset(1.3); 
   
   gPad->Modified();
   
   c1 -> SaveAs("btag_eff_83x_ttbar_hadron_semilep.png");
   
}
