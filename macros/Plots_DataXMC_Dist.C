int Plots_DataXMC_Dist()
{
   gStyle->SetOptStat(0);
   const char * trg = "L1_SingleJet100";
   
   TFile * fm = new TFile(Form("distributions_qcd_%s.root",trg),"old");
   TFile * fd = new TFile(Form("distributions_data_%s.root",trg),"old");
   
   
   // pt
   
   TH1F * hm_pt = (TH1F*) fm->Get("ht_pt");
   hm_pt -> SetName("hm_pt");
   hm_pt -> SetMarkerColor(kBlack);
   hm_pt -> SetLineColor(kBlack);
   TH1F * hd_pt = (TH1F*) fd->Get("ht_pt");
   hd_pt -> SetName("hd_pt");
   hd_pt -> SetMarkerColor(kRed);
   hd_pt -> SetLineColor(kRed);
   
   hm_pt -> Rebin(2);
   hd_pt -> Rebin(2);
     
   hm_pt -> Scale(1./hm_pt -> Integral());
   hd_pt -> Scale(1./hd_pt->GetEntries());
   
   std::cout << hm_pt -> Integral() << std::endl;
   std::cout << hd_pt -> Integral() << std::endl;
   
   hm_pt -> GetXaxis() -> SetRangeUser(20,200);
   hm_pt -> GetYaxis() -> SetTitle("entries (a.u.)");
   hm_pt -> GetXaxis() -> SetTitle("L1 leading jet pt (GeV)");
   hm_pt -> GetYaxis() -> SetTitleOffset(2);
   hm_pt -> SetTitle(Form("%s",trg)); 
    
   TCanvas * c_pt = new TCanvas("c_pt","",700,600);
   c_pt -> SetLogy();
   c_pt -> SetLeftMargin(0.15);
   
   hm_pt -> Draw("");
   hd_pt -> Draw("same");
   
   TLegend * leg_pt = new TLegend(0.65,0.70,0.85,0.80);   
   leg_pt->AddEntry(hm_pt,"QCD","p");
   leg_pt->AddEntry(hd_pt,"Data","p");
   
   
   leg_pt -> Draw();
   
   c_pt -> SaveAs(Form("%s_DataXMC_pt.png",trg));
   
   
   // eta
   
   TH1F * hm_eta = (TH1F*) fm->Get("ht_eta");
   hm_eta -> SetName("hm_eta");
   hm_eta -> SetMarkerColor(kBlack);
   hm_eta -> SetLineColor(kBlack);
   TH1F * hd_eta = (TH1F*) fd->Get("ht_eta");
   hd_eta -> SetName("hd_eta");
   hd_eta -> SetMarkerColor(kRed);
   hd_eta -> SetLineColor(kRed);
   
   hm_eta -> Scale(1./hm_eta -> Integral());
   hd_eta -> Scale(1./hd_eta->GetEntries());
   
   hm_eta  -> GetYaxis() -> SetRangeUser(0,0.05);
   hm_eta -> GetYaxis() -> SetTitle("entries (a.u.)");
   hm_eta -> GetXaxis() -> SetTitle("L1 leading jet eta");
   hm_eta -> GetYaxis() -> SetTitleOffset(2);
   hm_eta -> SetTitle(Form("%s",trg)); 
   
   
   TCanvas * c_eta = new TCanvas("c_eta","",700,600);
   c_eta -> SetLeftMargin(0.15);
   
   hm_eta -> Draw();
   hd_eta -> Draw("same");
   
   TLegend * leg_eta = new TLegend(0.65,0.70,0.85,0.80);   
   leg_eta->AddEntry(hm_eta,"QCD","p");
   leg_eta->AddEntry(hd_eta,"Data","p");
   
   
   leg_eta -> Draw();
   
   c_eta -> SaveAs(Form("%s_DataXMC_eta.png",trg));
   
   return 0;
}
