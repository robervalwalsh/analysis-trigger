int Efficiency()
{
   TFile * f_ref = new TFile("histograms_pfjet_jetht_2017E_hlt60_refhlt60_pt60to500_all.root","OLD");
   TFile * f_nom = new TFile("histograms_pfjet_jetht_2017E_hlt100_refhlt60_pt60to500_all.root","OLD");
   
   TH1F * h_ref = (TH1F*) f_ref -> Get("pt");
   TH1F * h_nom = (TH1F*) f_nom -> Get("pt");
   
   for ( int i = 0 ; i < 300 ; ++i )
   {
      if ( h_ref->GetBinContent(i+1) < h_nom->GetBinContent(i+1) )
      {
         h_ref->SetBinContent(i+1, h_nom->GetBinContent(i+1));
         h_ref->SetBinError(i+1, h_nom->GetBinError(i+1));
         std::cout << "bin " << i+1 << ": " <<  h_ref->GetBinContent(i+1) << "   " << h_nom->GetBinContent(i+1) << std::endl;
      }

   }
   
   
   
   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_nom,h_ref,"cl=0.683 b(1,1) mode");
//   TGraphAsymmErrors * g_eff = new TGraphAsymmErrors(h_ref,h_nom,"pois");
   
   TCanvas * c1 = new TCanvas("c1","");
   g_eff -> Draw("AP");
 
//   c1 -> SaveAs("eff_l1mu7.png");
   
   TFile * out = new TFile("eff_2017E.root","recreate");
   g_eff -> Write();
   
   return 0;
}
