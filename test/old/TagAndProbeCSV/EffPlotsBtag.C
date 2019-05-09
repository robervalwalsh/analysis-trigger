int EffPlotsBtag(std::string ptbin)
{
   std::string era[5] = {"2017C","2017D","2017E-v1","2017E-v2","2017F"};
   Color_t mycolor[5] = {kBlack,kRed,kBlue,kGreen,kMagenta};
   TFile * f[5];
   
   TCanvas * c1 = new TCanvas("c1","");
   TLegend * legend = new TLegend(0.6,0.2,0.8,0.45,"Efficiency per era");
   
   TMultiGraph *mg = new TMultiGraph();
   TGraphAsymmErrors * eff[5];
   for ( int i = 0; i < 5 ; ++i )
   {
      f[i] = new TFile(Form("eff_btag_%s_%s.root",era[i].c_str(),ptbin.c_str()),"OLD");
      eff[i] = (TGraphAsymmErrors*) f[i] -> Get("divide_btaglog_probe_num_cumulative_by_btaglog_probe_den_cumulative");
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
   mg -> GetYaxis()->SetNdivisions(510);
   mg -> GetXaxis()->SetNdivisions(510);
   mg -> GetXaxis()->SetRangeUser(0,10);
   mg -> GetYaxis()->SetRangeUser(0,1);
   mg -> GetXaxis()->SetTitle("-ln(1-CSVv2)");
   mg -> GetYaxis()->SetTitle("efficiency");
   for ( int i = 0; i < 5 ; ++i )
   {
      legend-> AddEntry(eff[i]->GetName(),era[i].c_str(),"ep");
   }
   legend -> Draw();
   
   c1->Modified();
//    c1->cd();
//    c1->SetSelected(c1);
//    c1->ToggleToolBar();
   
   c1 -> SaveAs(Form("OnlineBtag_EffBtag_%s_eras.png",ptbin.c_str()));
   
   
   return 0;
}
