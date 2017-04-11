void Plots_DataXDataPoint()
{
   int version = 0;
   
   int nbx=2544;
   const char * trg = "L1_DoubleJetC100Eta2p3_dEtaMax1p6";
//   trg = "L1_Mu10_DiJet32";
   
   // data point
   int n_dp = 4;
   double rate_dp[100];
   double ratee_dp[100];
   double pu_dp[100];
   double pue_dp[100];
   
   pu_dp[0] = 45;
   pu_dp[1] = 47;
   pu_dp[2] = 50;
   pu_dp[3] = 55;
   pue_dp[0] = 0;
   pue_dp[1] = 0;
   pue_dp[3] = 0;
   pue_dp[0] = 0;
   
   if ( std::string(trg) == "L1_Mu10_DiJet32" )
   {
      rate_dp [0] = 2254.95;
      ratee_dp[0] = 258.661;
      rate_dp [1] = 2366.78;
      ratee_dp[1] = 298.186;
      rate_dp [2] = 2669.37;
      ratee_dp[2] = 326.115;
      rate_dp [3] = 3698.15;
      ratee_dp[3] = 438.89;
   }
   if ( std::string(trg) == "L1_DoubleJetC100Eta2p3_dEtaMax1p6" )
   {
      rate_dp [0] = 3857.15;
      ratee_dp[0] = 338.295;
      rate_dp [1] = 4357.88;
      ratee_dp[1] = 404.619;
      rate_dp [2] = 4900.48;
      ratee_dp[2] = 441.861;
      rate_dp [3] = 4948.23;
      ratee_dp[3] = 507.678;
   }
   
   TGraphErrors * g1 = new TGraphErrors(n_dp,pu_dp,rate_dp,pue_dp,ratee_dp);
   g1 -> SetMarkerStyle(20);
   g1 -> SetMarkerColor(kBlue);
   g1 -> SetLineColor(kBlue);
   
   TFile * f2 = new TFile(Form("data_%d_bx_%s_offlumi.root",nbx,trg),"old");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   float minFit = 20;
   if ( std::string(trg) == "L1_DoubleJetC100" || std::string(trg) == "L1_DoubleJetC100Eta2p4" || std::string(trg) == "L1_Mu10_DiJet40" || std::string(trg) == "L1_Mu10_DiJet32" ) minFit = 0;
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);
   c1->SetLeftMargin(0.15);
   
   TF1 * pol2_2 = new TF1("pol2_2","pol2",0,70);
   pol2_2 -> SetLineColor(kRed);
   TFitResultPtr fit2 = g2 -> Fit(pol2_2,"S","",minFit,70);
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> SetTitle(Form("%s",trg));
   if ( std::string(trg) == "L1_Mu10_DiJet32" )
      mg -> SetTitle("L1_Mu10_dRMax0p4_DoubleJet32_dEtaMax1p6");  
   if ( std::string(trg) == "L1_Mu12_DiJet32" )
      mg -> SetTitle("L1_Mu12_dRMax0p4_DoubleJet32_dEtaMax1p6");  
   if ( std::string(trg) == "L1_Mu10_DiJet40" )
      mg -> SetTitle("L1_Mu10_dRMax0p4_DoubleJet40_dEtaMax1p6");  
   std::string gtitle = std::string(mg->GetTitle());
   if ( version == 1 )
   {
      mg -> SetTitle("CMS Work in Progress");
   }
   mg -> Add(g1);
   mg -> Add(g2);
   
   mg -> Draw("ap");
   
//    TH1F * h_c1 = c1->DrawFrame(20,0,45,1000);
//    h_c1 -> GetYaxis() -> SetTitleOffset(2);
//    h_c1 -> GetYaxis() -> SetTitle("rate (Hz)");
//    h_c1 -> GetXaxis() -> SetTitle("pile up");
//    
//    h_c1 -> Draw();
//    g1 -> Draw("ap");
//    g2 -> Draw("p");
   
   TLegend * leg = new TLegend(0.15,0.65,0.9,0.9); 
   leg->SetTextSize(0.031);
   if ( version == 0 )
   {
      leg->SetY1(0.8);
      leg->AddEntry(g1,"Data run 283171 (L1T Ntuple)","p");
      leg->AddEntry(g2,"Data parked ZB (pol2 fit)","lp");
   }
   
   
   leg -> Draw();
   
   mg -> GetYaxis() -> SetRangeUser(0,17E3);
   mg -> GetXaxis() -> SetRangeUser(20,62);
   mg -> GetYaxis() -> SetTitle("rate (Hz)");
   mg -> GetXaxis() -> SetTitle("pile up");
   mg -> GetYaxis() -> SetTitleOffset(2); 
     
//   h_c1 -> GetXaxis() -> SetRangeUser(20,65);
   
//   mg -> Draw("a");
//    mg -> SetMinimum(0);
//    mg -> SetMaximum(10000);
//    
   gPad->Modified();
   
  c1 -> SaveAs(Form("%s_DataXDataPoint_%d_bx_v%d.png",trg,nbx,version));  
 
}
