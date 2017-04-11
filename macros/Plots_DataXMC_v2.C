// data point
int n_dp = 4;
double rate_dp[100];
double ratee_dp[100];
double pu_dp[100];
double pue_dp[100];

void datapoints(const char *);

void Plots_DataXMC_v2()
{
   int version = 0;
   
   int nbx=2544;
   const char * trg = "L1_DoubleJetC100Eta2p4";
   datapoints(trg);
   
   
   TFile * f1 = new TFile(Form("qcd_%d_bx_%s.root",nbx,trg),"old");
   TFile * f2 = new TFile(Form("data_%d_bx_%s_offlumi.root",nbx,trg),"old");
   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g1 -> SetMarkerStyle(20);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   TGraphErrors * g3 = new TGraphErrors(n_dp,pu_dp,rate_dp,pue_dp,ratee_dp);
   g3 -> SetMarkerStyle(20);
   g3 -> SetMarkerColor(kBlue);
   g3 -> SetLineColor(kBlue);
   
   
   
   float minFit = 20;
   if ( std::string(trg) == "L1_DoubleJetC100" || std::string(trg) == "L1_DoubleJetC100Eta2p4" || std::string(trg) == "L1_Mu10_DiJet40" || std::string(trg) == "L1_Mu10_DiJet32" ) minFit = 0;
   
   TCanvas * c1 = new TCanvas("c1","", 750,600);
   c1->SetLeftMargin(0.15);
   
   TF1 * pol2_1 = new TF1("pol2_1","pol2",0,70);
   pol2_1 -> SetLineColor(kBlack);
   TFitResultPtr fit1 = g1 -> Fit(pol2_1,"S","",minFit,70);
   
   TF1 * pol2_2 = new TF1("pol2_2","pol2",0,70);
   pol2_2 -> SetLineColor(kRed);
   TFitResultPtr fit2 = g2 -> Fit(pol2_2,"S","",minFit,70);
   
   TMultiGraph * mg = new TMultiGraph();
   
   mg -> SetTitle(Form("%s",trg));
   if ( std::string(trg) == "L1_Mu10_DiJet32" )
      mg -> SetTitle("L1_Mu10_Jet32Eta2p3_dRMax0p4_DoubleJet32Eta2p3_dEtaMax1p6");  
   if ( std::string(trg) == "L1_Mu12_DiJet32" )
      mg -> SetTitle("L1_Mu12_dRMax0p4_DoubleJet32_dEtaMax1p6");  
   if ( std::string(trg) == "L1_Mu10_DiJet40" )
      mg -> SetTitle("L1_Mu10_dRMax0p4_DoubleJet40_dEtaMax1p6");  
   if ( std::string(trg) == "L1_DoubleJetC100Eta2p4" )
      mg -> SetTitle("L1_DoubleJet100Eta2p4");  
   std::string gtitle = std::string(mg->GetTitle());
   if ( version == 1 )
   {
      mg -> SetTitle("CMS Work in Progress");
   }
   mg -> SetTitle("");
   mg -> Add(g1);
   mg -> Add(g2);
   mg -> Add(g3);
   
   mg -> Draw("ap");
   
//    TH1F * h_c1 = c1->DrawFrame(20,0,45,1000);
//    h_c1 -> GetYaxis() -> SetTitleOffset(2);
//    h_c1 -> GetYaxis() -> SetTitle("rate (Hz)");
//    h_c1 -> GetXaxis() -> SetTitle("pile up");
//    
//    h_c1 -> Draw();
//    g1 -> Draw("ap");
//    g2 -> Draw("p");
   
   TLegend * leg = new TLegend(0.15,0.75,0.9,0.9); 
   leg->SetTextSize(0.031);
   if ( version == 0 )
   {
      leg->SetY1(0.70);
      leg->SetHeader(gtitle.c_str());
      leg->AddEntry(g1,"80X TSG QCD Monte Carlo (pol2 fit)","lp");
      leg->AddEntry(g2,"2016 data parked ZB (pol2 fit)","lp");
      leg->AddEntry(g3,"2016 data run 283171 L1T Ntuple","p");
   }
   if ( version == 1 )
   {
      leg->SetHeader(gtitle.c_str());
      leg->AddEntry(g2,"pp 13 TeV collisions data","p");
      leg->AddEntry(pol2_2,"Pol2 fit to pp 13 TeV collisions data","l");
      leg->AddEntry(g1,"QCD Monte Carlo pp 13 TeV","p");
      leg->AddEntry(pol2_1,"Pol2 fit to QCD Monte Carlo pp 13 TeV","l");
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
   
  c1 -> SaveAs(Form("%s_DataXMC_v2_%d_bx_v%d.png",trg,nbx,version));  
 
}


void datapoints(const char * trg )
{
   
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
   
   if ( std::string(trg) == "L1_DoubleJetC100" )
   {
      rate_dp [0] = 6794.52;
      ratee_dp[0] = 448.995;
      rate_dp [1] = 6536.82;
      ratee_dp[1] = 495.555;
      rate_dp [2] = 8366.67;
      ratee_dp[2] = 577.355;
      rate_dp [3] = 8698.48;
      ratee_dp[3] = 673.108;
   }
   
   if ( std::string(trg) == "L1_DoubleJetC112" )
   {
      rate_dp [0] = 4302.21;
      ratee_dp[0] = 357.279;
      rate_dp [1] = 3944.63;
      ratee_dp[1] = 384.957;
      rate_dp [2] = 5059.84;
      ratee_dp[2] = 448.988;
      rate_dp [3] = 5312.84;
      ratee_dp[3] = 526.05;
   }
   
   if ( std::string(trg) == "L1_DoubleJetC100Eta2p4" )
   {
      rate_dp [0] = 5904.41;
      ratee_dp[0] = 418.553;
      rate_dp [1] = 6161.14;
      ratee_dp[1] = 481.104;
      rate_dp [2] = 7569.84;
      ratee_dp[2] = 549.174;
      rate_dp [3] = 7708.83;
      ratee_dp[3] = 633.662;
   }
   
}
