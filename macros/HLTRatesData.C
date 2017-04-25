void HLTRatesData(const std::string & trigger, const std::string & file)
{
   // All hadronic
//   trigger = "HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw";
//   trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw";
//   trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_psw";
//   trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5_psw";
   
//   trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw";
//   trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw";
//   trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_psw";
//   trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5_psw";
   
   

   double nbx = 2600;
   TFile * inFile = new TFile(file.c_str(),"old");
//   TFile * inFile = new TFile("../test/80x/BTagCSV/mssmhbb_triggers_data_BTagCSV.root","old");
   TH1F * hLumis   = (TH1F*) inFile -> Get("Lumis");
   TH1F * hpath = (TH1F*) inFile -> Get(Form("h_n%s",trigger.c_str()));
   
   double pileup[100];
   double pileupErr[100];
   double rate[100];
   double rateErr[100];
   int n = 0;
   
   cout << hLumis->GetNbinsX() << endl;
   for ( int i = 0; i < hLumis->GetNbinsX(); ++i )
   {
      double pu = hLumis->GetBinCenter(i+1);
      if ( pu < 14 || pu == 32.5 || pu > 44 ) continue;
      pileup[n] = pu;
      pileupErr[n] = 0;
      double nlumi = hLumis->GetBinContent(i+1);
      double nevts = hpath->GetBinContent(i+1);
      double nevtsErr = hpath->GetBinError(i+1);
      rate[n] = nevts*nbx/nlumi/23.31;
      rateErr[n]= nevtsErr*nbx/nlumi/23.31;
//      if ( nlumi == 0 || pu > 42 )  break;
      
      ++n;
      
      
   }
   
   TGraphErrors * gRate = new TGraphErrors(n,pileup,rate,pileupErr,rateErr);
   gRate -> SetName("rate_total");
   gRate -> SetMarkerStyle(20);
   gRate -> Draw("AP");
   
   TFile * out = new TFile(Form("data_%d_bx_%s.root",int(nbx),trigger.c_str()),"recreate");
   gRate -> Write();
   
   out -> Close();

   
}

