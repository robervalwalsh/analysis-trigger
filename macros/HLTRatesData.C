void HLTRatesData(const std::string & trigger, const std::string & file)
{
   // All hadronic
//   trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw";
   
//   trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw";
   
// root -l HLTRatesData.C\(\"HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw\",\"../test/80x/BTagMu20/mssmhbb_triggers_data_BTagMu20.root\"\)   
// root -l HLTRatesData.C\(\"HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw\",\"../test/80x/BTagCSV/mssmhbb_triggers_data_BTagCSV.root\"\) 

//root -l HLTRatesData.C\(\"HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350_psw\",\"../test/80x/BTagMu20/mssmhbb_triggers_data_BTagMu20_control_backup_v2.root\"\)
//root -l HLTRatesData.C\(\"HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130_psw\",\"../test/80x/BTagCSV/mssmhbb_triggers_data_BTagCSV_control_backup_v2.root\"\)
//root -l HLTRatesData.C\(\"HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350_psw\",\"../test/80x/BTagCSVControl/mssmhbb_triggers_data_BTagCSVControl_control_backup_v2.root\"\)

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
      if ( pu < 14 || pu > 43 ) continue;
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
   
   TFile * out = new TFile(Form("data2016_%d_bx_%s.root",int(nbx),trigger.c_str()),"recreate");
   gRate -> Write();
   
   out -> Close();

   
}

