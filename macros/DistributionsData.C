int DistributionsData()
{
   const char * trg = "L1_DoubleJetC112";
   
   TH1::SetDefaultSumw2();
   std::vector<TFile *> f;
   
   f.push_back(new TFile("../test/80x/mssmhbb_triggers_distributions_data_parking_zerobias0_offlumi.root","old"));
   
   TH1F * ht_n   ;
   TH1F * ht_pt  ;
   TH1F * ht_eta ;
   TH1F * ht_phi ;
   
   for ( size_t i = 0 ; i < f.size() ; ++i )  // loop over files
   {
      TH1F * h_n    = (TH1F*) f[i]->Get(Form("h_jetsN_%s",trg));
      TH1F * h_pt   = (TH1F*) f[i]->Get(Form("h_jet1PT_%s",trg));
      TH1F * h_eta  = (TH1F*) f[i]->Get(Form("h_jet1ETA_%s",trg));
      TH1F * h_phi  = (TH1F*) f[i]->Get(Form("h_jet1PHI_%s",trg));
      
      if ( i == 0 )
      {
         ht_n    =  (TH1F*)   h_n     -> Clone("ht_n");
         ht_pt   =  (TH1F*)   h_pt    -> Clone("ht_pt");
         ht_eta  =  (TH1F*)   h_eta   -> Clone("ht_eta");
         ht_phi  =  (TH1F*)   h_phi   -> Clone("ht_phi");
         ht_n    -> Reset();
         ht_pt   -> Reset();
         ht_eta  -> Reset();
         ht_phi  -> Reset();
      }
      
      ht_n    ->  Add(h_n  );
      ht_pt   ->  Add(h_pt );
      ht_eta  ->  Add(h_eta);
      ht_phi  ->  Add(h_phi);
      
   }
   
   ht_pt  -> SetMarkerStyle(20);
   ht_eta -> SetMarkerStyle(20);
   ht_phi -> SetMarkerStyle(20);
   
   TFile * fout = new TFile(Form("distributions_data_%s.root",trg),"recreate");
   
   ht_pt -> Write();
   ht_eta -> Write();
   ht_phi -> Write();
   
   
   return 0;
   
}
