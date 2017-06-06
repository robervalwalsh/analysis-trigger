int plot()
{

   TFile * f = new TFile("mssmhbb_triggers_data_test.root","old");
   
   TH1F * h[4];
   
   h[0]  = (TH1F*) f -> Get("muon12");
   h[1]  = (TH1F*) f -> Get("l1muon12");
   h[2]  = (TH1F*) f -> Get("l2muon8");
   h[3]  = (TH1F*) f -> Get("l3muon12");

   int rebin = 2;
   
   TCanvas * c1 = new TCanvas("c1","");
   
   c1 -> SetLogy();
   
   int color[4] = {kBlack,kGreen,kRed,kBlue};
   
   h[0] -> GetXaxis() -> SetTitle("muon pt (GeV)");
   h[0] -> GetYaxis() -> SetTitle("entries");
   
   
   for ( int i = 0; i < 4 ; ++i )
   {
       h[i] -> Rebin(rebin);
       h[i] -> SetMarkerStyle(20);
       h[i] -> SetMarkerColor(color[i]);
       
       if ( i == 0 ) h[i] -> Draw();
       else if ( i > 1 ) h[i] -> Draw("same");
   }
   
   c1 -> SaveAs("muon_pt.png");
   
   return 0;

}
