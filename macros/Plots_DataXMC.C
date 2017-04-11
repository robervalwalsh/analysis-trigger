void Plots_DataXMC()
{
   int version = 0;
   
   int nbx=2544;
   const char * trg = "L1_Mu10_DiJet32";
   
   TFile * f1 = new TFile(Form("qcd_%d_bx_%s.root",nbx,trg),"old");
   TFile * f2 = new TFile(Form("data_%d_bx_%s_offlumi.root",nbx,trg),"old");
   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g1 -> SetMarkerStyle(20);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   float minFit = 20;
   if ( std::string(trg) == "L1_DoubleJetC100" || std::string(trg) == "L1_DoubleJetC100Eta2p4" || std::string(trg) == "L1_Mu10_DiJet40" || std::string(trg) == "L1_Mu10_DiJet32" ) minFit = 0;
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);
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
      leg->AddEntry(g1,"QCD (pol2 fit)","lp");
      leg->AddEntry(g2,"Data (pol2 fit)","lp");
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
   
  c1 -> SaveAs(Form("%s_DataXMC_%d_bx_v%d.png",trg,nbx,version));  
 
}
