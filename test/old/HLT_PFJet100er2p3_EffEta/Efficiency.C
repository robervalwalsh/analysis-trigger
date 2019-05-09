int Efficiency()
{
   std::vector<std::string> erasList = {"2017C","2017D","2017E","2017F","2017C-v1","2017C-v2","2017C-v3","2017E-v1","2017E-v2"};
   std::vector<Color_t> colorsList = {kBlack,kRed,kBlue,kMagenta,kBlack,kBlack+2,kGray,kBlue,kCyan};
   std::map<std::string,Color_t> colors;
   std::map<std::string, std::string> files;
   for ( auto & era : erasList )
   {
      auto i = &era - &erasList[0];
      colors[era] = colorsList[i];
      std::string myfile = Form("histograms_pfjet_jetht_%s_HLT_PFJet100er2p3_refHLT_PFJet100er3p0_all.root",era.c_str());
      if (std::ifstream(myfile))
      {
         files[era] = myfile;
      }
   }
   
   TCanvas * c1 = new TCanvas("c1","");
   TLegend * legend = new TLegend(0.65,0.2,0.85,0.45,"Efficiency per era");
   
   std::map< std::string, TGraphAsymmErrors *> eff;
   TMultiGraph *mg = new TMultiGraph();
   
   for ( auto & f : files )
   {
      std::string era = f.first;
      TFile * tf = new TFile(f.second.c_str(),"old");
      TH1F * h_den = (TH1F*) tf -> Get("eta_den");
      TH1F * h_num = (TH1F*) tf -> Get("eta_num");
      eff[era] = new TGraphAsymmErrors(h_num,h_den,"cl=0.683 b(1,1) mode");
      eff[era] -> SetMarkerStyle(20);
      eff[era] -> SetMarkerSize(0.8);
      eff[era] -> SetMarkerColor(colors[era]);
      eff[era] -> SetLineColor(colors[era]);
      eff[era] -> SetName(Form("eff_%s",era.c_str()));
      mg -> Add(eff[era]);
   }
   
   c1->SetGridx();
   c1->SetGridy();
   c1->SetTickx(1);
   c1->SetTicky(1);

   mg -> Draw("ap");
   mg -> GetYaxis()->SetNdivisions(520);
   mg -> GetXaxis()->SetNdivisions(520);
   mg -> GetXaxis()->SetRangeUser(-2.5,2.5);
   mg -> GetXaxis()->SetTitle("jet eta");
   mg -> GetYaxis()->SetTitle("efficiency");

   for ( auto & e : eff )
      legend-> AddEntry(e.second->GetName(),e.first.c_str(),"ep");

   legend -> Draw();
   
   c1->Modified();
   
   c1 -> SaveAs("PFJet100_EffEta_eras.png");
   

   return 0;
}
