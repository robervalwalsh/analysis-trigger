int sf_fit()
{
    bool compare = true;
    float xmin = 40;
    float xmax = 1000;
    float ymin = 0.3;
    float ymax = 1.7;
    int nbins = int(xmax - xmin);


    TFile * f = new TFile("online_btag_SF_j3veto_probbb_2018.root","old");
    TGraphAsymmErrors * g_sf = (TGraphAsymmErrors*) f->Get("scale_factors_graph");
    g_sf -> SetMarkerColor(kRed);
    g_sf -> SetLineColor(kRed);

    g_sf->GetFunction("fpol")->Delete();

    // g_sf -> RemovePoint(g_sf->GetN()-1);

    TCanvas * c1 = new TCanvas("c1","");

    // char const * const x_shift = "0";
    TF1 *fpol = new TF1("fpol", "pol1", xmin, xmax);
    fpol -> SetLineColor(kRed);
    fpol -> SetLineStyle(2);
    auto fit_result = g_sf->Fit(fpol, "QSR");
    fpol -> FixParameter(0,fpol->GetParameter(0));
    // fpol -> FixParameter(1,fpol->GetParameter(1));
    fit_result = g_sf->Fit(fpol, "QSR");
    fit_result ->Print("V");


    //Create a histo with errors to store the intervals
    int ngr = g_sf->GetN();
    TH1D *hint = new TH1D("hint","", nbins, xmin, xmax);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,0.68);
    hint->SetFillColor(kRed);
    hint->SetFillStyle(3002);
    hint->GetYaxis()->SetRangeUser(ymin,ymax);
    hint->GetYaxis()->SetTitle("rel. online btag scale factor");
    hint->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    hint->SetStats(false);
    // hint->GetXaxis()->SetRangeUser(xmin,xmax);
    hint ->SetLineColor(kRed);
    hint->Draw("e3");

    g_sf->Draw("psame");


    TGraphAsymmErrors * g_sf_2 = nullptr;
    TH1D *hint_2 = nullptr;
    TF1 *fpol_2 = nullptr;
    TFitResultPtr fit_result_2 = nullptr;
    if ( compare )
    {
        TFile * f_2 = new TFile("online_btag_SF__2018.root","old");
        g_sf_2 = (TGraphAsymmErrors*) f_2->Get("scale_factors_graph");
        g_sf_2 -> GetFunction("fpol") ->Delete();
        // g_sf_2 -> RemovePoint(g_sf_2->GetN()-1);

        g_sf_2 -> SetMarkerStyle(20);
        g_sf_2 -> SetMarkerColor(kBlue);
        g_sf_2 -> SetLineColor(kBlue);
        fpol_2 = new TF1("fpol_2", "pol1", xmin, xmax);
        fpol_2 -> SetLineColor(kBlue);
        fpol_2 -> SetLineStyle(2);
        fit_result_2 = g_sf_2->Fit(fpol_2, "QSR");
        fpol_2 -> FixParameter(0,fpol_2->GetParameter(0));
        // fpol_2 -> FixParameter(1,fpol_2->GetParameter(1));
        fit_result_2 = g_sf_2->Fit(fpol_2, "QSR");
        //Create a histo with errors to store the intervals
        int ngr_2 = g_sf_2->GetN();
        hint_2 = new TH1D("hint_2","", nbins, xmin, xmax);
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint_2,0.68);
        hint_2->SetFillColor(kBlue);
        hint_2->SetFillStyle(3002);
        hint_2->GetYaxis()->SetRangeUser(ymin,ymax);
        hint_2->GetYaxis()->SetTitle("rel. online btag scale factor");
        hint_2->GetXaxis()->SetTitle("jet p_{T} [GeV]");
        hint_2->SetStats(false);
        hint_2->GetXaxis()->SetRangeUser(xmin,xmax);
        hint_2->Draw("e3 same");
        g_sf_2->Draw("psame");
        fit_result_2 ->Print("V");
    }

    c1 -> SaveAs("online_btag_SF_j3veto_probbb_2018_compare_fixpars_linear.png");
    return 0;
}