int EfficiencyEta()
{
   TFile * f = new TFile("histograms_2017F_tep_400to1000_all.root","OLD");
   
   TH1F * h_ref = (TH1F*) f -> Get("eta_probe_den");
   TH1F * h_nom = (TH1F*) f -> Get("eta_probe_num");
   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_nom,h_ref,"cl=0.683 b(1,1) mode");
//   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_ref,h_nom,"pois");
   
   TCanvas * c1 = new TCanvas("c1","");
   g_eff -> Draw("AP");
 
//   c1 -> SaveAs("eff_l1mu7.png");
   
   TFile * out = new TFile("eff_eta_2017F_pt_400to1000.root","recreate");
   g_eff -> Write();
   
   return 0;
}
