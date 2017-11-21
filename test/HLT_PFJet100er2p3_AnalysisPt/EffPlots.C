int EffPlots()
{
   std::string era[4] = {"2017C","2017D","2017E","2017F"};
   Color_t mycolor[4] = {kBlack,kRed,kBlue,kMagenta};
   TFile * f[4];
   
   TCanvas * c1 = new TCanvas("c1","");
   TLegend * legend = new TLegend(0.65,0.2,0.85,0.45,"Efficiency per era");
   
   TMultiGraph *mg = new TMultiGraph();
   TGraphAsymmErrors * eff[4];
   for ( int i = 0; i < 4 ; ++i )
   {
      f[i] = new TFile(Form("eff_%s.root",era[i].c_str()),"OLD");
      eff[i] = (TGraphAsymmErrors*) f[i] -> Get("divide_pt_by_pt");
      eff[i] -> SetMarkerStyle(20);
      eff[i] -> SetMarkerSize(0.8);
      eff[i] -> SetMarkerColor(mycolor[i]);
      eff[i] -> SetLineColor(mycolor[i]);
      eff[i] -> SetName(Form("eff_%d",i));
      mg -> Add(eff[i]);
   }
   
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);

   mg -> Draw("ap");
   mg -> GetYaxis()->SetNdivisions(520);
   mg -> GetXaxis()->SetNdivisions(510);
   mg -> GetXaxis()->SetRangeUser(55,250);
   mg -> GetXaxis()->SetTitle("jet PT (GeV)");
   mg -> GetYaxis()->SetTitle("efficiency");
   for ( int i = 0; i < 4 ; ++i )
   {
      legend-> AddEntry(eff[i]->GetName(),era[i].c_str(),"ep");
   }
   legend -> Draw();
   
   c1->Modified();
//    c1->cd();
//    c1->SetSelected(c1);
//    c1->ToggleToolBar();
   
   c1 -> SaveAs("PFJet100_EffPt_eras.png");
   
   
   return 0;
}
