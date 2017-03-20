void PlotsData_Lumi_OnXOff()
{
   TFile * f1 = new TFile("data_2208-bx_L1_SingleJet35_onlumi.root","old");
   TFile * f2 = new TFile("data_2208-bx_L1_SingleJet35_offlumi.root","old");
   TGraphAsymmErrors * g1 = (TGraphAsymmErrors*) f1 -> Get("rate_total");
   TGraphAsymmErrors * g2 = (TGraphAsymmErrors*) f2 -> Get("rate_total");
   g1 -> SetMarkerStyle(20);
   g2 -> SetMarkerStyle(20);
   g2 -> SetMarkerColor(kRed);
   g2 -> SetLineColor(kRed);
   
   TCanvas * c1 = new TCanvas("c1","", 700,600);
   c1->SetLeftMargin(0.15);
   
   TF1 * pol2_1 = new TF1("pol2_1","pol2",20,70);
   pol2_1 -> SetLineColor(kBlack);
   TFitResultPtr fit1 = g1 -> Fit(pol2_1,"S","",20,70);
   
   TF1 * pol2_2 = new TF1("pol2_2","pol2",20,70);
   pol2_2 -> SetLineColor(kRed);
   TFitResultPtr fit2 = g2 -> Fit(pol2_2,"S","",20,70);
   
   TMultiGraph * mg = new TMultiGraph();
   
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
   
   TLegend * leg = new TLegend(0.15,0.8,0.9,0.9);   
   leg->AddEntry(g1,"L1_SingleJet35 - online lumi (pol2 fit)","lp");
   leg->AddEntry(g2,"L1_SingleJet35 - offline lumi (pol2 fit)","lp");
   
   
   leg -> Draw();
   
   mg -> GetYaxis() -> SetRangeUser(0,6E6);
   mg -> GetXaxis() -> SetRangeUser(20,45);
   mg -> GetYaxis() -> SetTitle("rate (Hz)");
   mg -> GetXaxis() -> SetTitle("pile up");
   mg -> GetYaxis() -> SetTitleOffset(2); 
     
//   h_c1 -> GetXaxis() -> SetRangeUser(20,65);
   
//   mg -> Draw("a");
//    mg -> SetMinimum(0);
//    mg -> SetMaximum(10000);
//    
   gPad->Modified();
   
  c1 -> SaveAs("L1_SingleJet35_Data_Lumi_OnXOff.png");  
 
}
