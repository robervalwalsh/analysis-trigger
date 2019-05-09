int Efficiency(std::string mydir, std::string era)
{
   TFile * f = new TFile(Form("%s/histograms_tep_%s.root",mydir.c_str(),era.c_str()),"OLD");
   
   TH1F * h_ref = (TH1F*) f -> Get("pt_probe_den_var");
   TH1F * h_nom = (TH1F*) f -> Get("pt_probe_num_var");
   
   h_ref->Rebin(1);
   h_nom->Rebin(1);
   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_nom,h_ref,"cl=0.683 b(1,1) mode");
//   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_ref,h_nom,"pois");
   
   TCanvas * c1 = new TCanvas("c1","");
   g_eff -> Draw("AP");
 
//   c1 -> SaveAs("eff_l1mu7.png");
   
   TFile * out = new TFile(Form("%s/eff_pt_%s.root",mydir.c_str(),era.c_str()),"recreate");
   g_eff -> Write();
   
   return 0;
}
