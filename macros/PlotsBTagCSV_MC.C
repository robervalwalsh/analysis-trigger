void PlotsBTagCSV_MC()
{
   TFile * f1 = new TFile("qcd_2600_bx_HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100.root","old");
   TFile * f2 = new TFile("qcd_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100.root","old");
   TFile * f3 = new TFile("qcd_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1.root","old");
   TFile * f4 = new TFile("qcd_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5.root","old");

   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   TGraphAsymmErrors * g3 = (TGraphAsymmErrors*) f3 -> Get("rate_total");
   TGraphAsymmErrors * g4 = (TGraphAsymmErrors*) f4 -> Get("rate_total");
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> Add(g1);
   mg -> Add(g2);
   mg -> Add(g3);
   mg -> Add(g4);
   
   mg -> Draw("ap");
   
}
