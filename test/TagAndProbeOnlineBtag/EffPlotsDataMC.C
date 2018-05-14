int EffPlotsDataMC(std::string mydir = ".")
{
//   std::string era[5] = {"2017C","2017D","2017E-v1","2017E-v2","2017F"};
//   Color_t mycolor[5] = {kBlack,kRed,kBlue,kGreen,kMagenta};
//   std::string era[5] = {"2017C","2017D","2017E","2017F","2017F"};
//   Color_t mycolor[5] = {kBlack,kRed,kBlue,kMagenta,kGreen};
   std::string era[5] = {"2017CDE","2017F","QCD","2017F","2017F"};
   std::string dirs[5];
   dirs[0] = "oldbin/ReReco17Nov2017_2017/deepcsv_0.8001_0.4941";
   dirs[1] = "oldbin/ReReco17Nov2017_2017/deepcsv_0.8001_0.4941";
   dirs[2] = "MonteCarlo/deepcsv_0.8001_0.4941";
   Color_t mycolor[5] = {kBlack,kRed,kBlue,kMagenta,kGreen};
   TFile * f[5];
   
   int neras = 3;
   
   TCanvas * c1 = new TCanvas("c1","");
   // top
   TLegend * legend = new TLegend(0.65,0.6,0.85,0.85,"Efficiency per era");
   // bottom
//   TLegend * legend = new TLegend(0.65,0.2,0.85,0.45,"Efficiency per era");
   
   TMultiGraph *mg = new TMultiGraph();
   TGraphAsymmErrors * eff[5];
   for ( int i = 0; i < neras ; ++i )
   {
      f[i] = new TFile(Form("%s/eff_pt_%s.root",dirs[i].c_str(),era[i].c_str()),"OLD");
      eff[i] = (TGraphAsymmErrors*) f[i] -> Get("divide_pt_probe_num_var_by_pt_probe_den_var");
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
   mg -> GetXaxis()->SetRangeUser(40,1000);
   mg -> GetYaxis()->SetRangeUser(0,1);
   mg -> GetXaxis()->SetTitle("jet PT (GeV)");
   mg -> GetYaxis()->SetTitle("efficiency");
   for ( int i = 0; i < neras ; ++i )
   {
      legend-> AddEntry(eff[i]->GetName(),era[i].c_str(),"ep");
   }
   legend -> Draw();
   
   c1->Modified();
//    c1->cd();
//    c1->SetSelected(c1);
//    c1->ToggleToolBar();
   
   c1 -> SaveAs("OnlineBtag_EffPt_DataMC.png");
   
   
   return 0;
}
