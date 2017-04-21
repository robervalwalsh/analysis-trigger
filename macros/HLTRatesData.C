void HLTRatesData()
{
   double nbx = 2500;
   TFile * inFile = new TFile("../test/80x/BTagMu20/mssmhbb_triggers_data_btagmu20.root","old");
//   TFile * inFile = new TFile("../test/80x/HighPU/mssmhbb_triggers_data_zerobiasbt.root","old");
   TH1F * hLumis   = (TH1F*) inFile -> Get("Lumis");
   TH1F * hsemilep = (TH1F*) inFile -> Get("h_nHLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5_psw");
   
   double pileup[100];
   double pileupErr[100];
   double rate[100];
   double rateErr[100];
   int n = 0;
   
   cout << hLumis->GetNbinsX() << endl;
   for ( int i = 0; i < hLumis->GetNbinsX(); ++i )
   {
      double pu = hLumis->GetBinCenter(i+1);
      if ( pu < 14 ) continue;
      pileup[n] = pu;
      pileupErr[n] = 0;
      double nlumi = hLumis->GetBinContent(i+1);
      double nevts = hsemilep->GetBinContent(i+1);
      double nevtsErr = hsemilep->GetBinError(i+1);
      rate[n] = nevts*nbx/nlumi/23.31;
      rateErr[n]= nevtsErr*nbx/nlumi/23.31;
//      if ( nlumi == 0 || pu > 42 )  break;
      
      ++n;
      
      
   }
   
   TGraphErrors * gRate = new TGraphErrors(n,pileup,rate,pileupErr,rateErr);
   gRate -> SetMarkerStyle(20);
   gRate -> Draw("AP");
   

}

