void PlotsBTagCSV()
{
   TFile * f1 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw.root","old");
   TFile * f2 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw.root","old");
   TFile * f3 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_psw.root","old");
   TFile * f4 = new TFile("data_2600_bx_HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5_psw.root","old");

   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("Graph");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("Graph");
   TGraphAsymmErrors * g3 = (TGraphAsymmErrors*) f3 -> Get("Graph");
   TGraphAsymmErrors * g4 = (TGraphAsymmErrors*) f4 -> Get("Graph");
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> Add(g1);
   mg -> Add(g2);
   mg -> Add(g3);
   mg -> Add(g4);
   
   mg -> Draw("ap");
   
}
