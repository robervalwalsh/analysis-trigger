int EfficiencyBtag(std::string era, std::string ptbin)
{
   TFile * f = new TFile(Form("histograms_%s_tep_%s_all.root",era.c_str(),ptbin.c_str()),"OLD");
   
   TH1F * h_ref = (TH1F*) f -> Get("btaglog_probe_den");
   TH1F * h_nom = (TH1F*) f -> Get("btaglog_probe_num");
   
//    TH1F* h_ref_cum  = (TH1F*) h_ref -> Clone("h_ref_cum");
//    h_ref_cum -> Reset();
//    TH1F* h_nom_cum  = (TH1F*) h_nom -> Clone("h_nom_cum");
//    h_nom_cum -> Reset();
   
   TH1F* h_ref_cum  = (TH1F*) h_ref -> GetCumulative(false);
   TH1F* h_nom_cum  = (TH1F*) h_nom -> GetCumulative(false);
   
   for ( int i = 0; i < h_ref_cum->GetNbinsX() ; ++i )
   {
//        if ( i < 1 )
//        {
//            h_ref_cum -> SetBinContent(i+1,0.);
//            h_nom_cum -> SetBinContent(i+1,0.);
//        }
       h_ref_cum -> SetBinError(i+1,sqrt(h_ref_cum->GetBinContent(i+1)));
       h_nom_cum -> SetBinError(i+1,sqrt(h_nom_cum->GetBinContent(i+1)));
   }
   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_nom_cum,h_ref_cum,"cl=0.683 b(1,1) mode");
//   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_ref,h_nom,"pois");
   
   TCanvas * c1 = new TCanvas("c1","");
   g_eff -> Draw("AP");
 
//   c1 -> SaveAs("eff_l1mu7.png");
   
   TFile * out = new TFile(Form("eff_btag_%s_%s.root",era.c_str(),ptbin.c_str()),"recreate");
   g_eff -> Write();
   
   return 0;
}
