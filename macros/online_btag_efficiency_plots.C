
int pt_selection(TH1F * h_pt, const float & min_pt, const float & max_pt);

int online_btag_efficiency_plots(const std::string & input_dir, const int & atype = 0)
{
    gStyle->SetOptTitle(0);
    std::vector<int> thresholds = {40,100,200,350};
    std::vector<int> pt_min = {30,100,200,350};
    std::vector<int> pt_max = {100,200,350,600};

    std::vector<int> colors = {kBlack, kRed, kBlue, kMagenta};
    std::map<int, TFile *> input_files;
    std::map<int, TH1F *> h_pt_probes;
    std::map<int, TH1F *> h_pt_probes_match;
    std::map<int, TH1F *> h_btag_jet1;
    std::map<int, TH1F *> h_btag_jet2;
    std::map<int, TGraphAsymmErrors *> g_efficiencies;

    TCanvas * c_efficiency = new TCanvas("c_efficiency","");
    TMultiGraph * mg_efficiency = new TMultiGraph();

    std::string type = "";
    if ( atype == 1 ) type = "_Mu12";

    for ( size_t i = 0; i <  thresholds.size() ; ++i )
    {
        auto thr = thresholds[i];
        input_files[thr] = new TFile(Form("%s/btag_probe_HLT%s_DoublePFJets%d_PNetBTag_0p11.root", input_dir.c_str(), type.c_str(), thr), "old");
        h_pt_probes[thr] = (TH1F*) input_files[thr] ->Get("pt_probe");
        h_pt_probes_match[thr] = (TH1F*) input_files[thr] ->Get("pt_probe_match");
        pt_selection(h_pt_probes[thr],pt_min[i],pt_max[i]);
        pt_selection(h_pt_probes_match[thr],pt_min[i],pt_max[i]);
        g_efficiencies[thr] = new TGraphAsymmErrors(h_pt_probes_match[thr],h_pt_probes[thr],"cl=0.683 b(1,1) mode");
        auto color = colors[i];
        g_efficiencies[thr] -> SetLineColor(color);
        g_efficiencies[thr] -> SetMarkerColor(color);
        g_efficiencies[thr] -> SetTitle(Form("HLT%s_DoublePFJets%d_PNetBTag_0p11",type.c_str(),thr));
        mg_efficiency->Add(g_efficiencies[thr]);
        // btag
        h_btag_jet1[thr] = (TH1F*) input_files[thr] ->Get("discr_offline_btag_jet1");
        h_btag_jet1[thr] -> SetLineColor(color);
        h_btag_jet2[thr] = (TH1F*) input_files[thr] ->Get("discr_offline_btag_jet2");
        h_btag_jet2[thr] -> SetLineColor(color);
    }
    mg_efficiency -> GetYaxis() -> SetTitle("relative online btag efficiency");
    mg_efficiency -> GetYaxis() -> SetRangeUser(0,1);
    mg_efficiency -> GetXaxis() -> SetTitle("probe jet pt (GeV)");

    c_efficiency ->cd();
    c_efficiency -> SetGrid();
    mg_efficiency ->Draw("ap");
    c_efficiency -> BuildLegend(0.40,0.55,0.88,0.88);

    std::string output_efficiency       = input_dir + "/" + input_dir + "_relative_online_btag_efficiency_fh.png";
    if ( atype == 1 ) output_efficiency = input_dir + "/" + input_dir + "_relative_online_btag_efficiency_sl.png";
    c_efficiency -> SaveAs(output_efficiency.c_str());

    TCanvas * c_btag_jet2 = new TCanvas("c_btag_jet2","");
    c_btag_jet2 -> SetGrid();
    for ( size_t i = 0; i <  thresholds.size() ; ++i )
    {
        auto thr = thresholds[i];
        h_btag_jet2[thr] -> GetYaxis() -> SetTitle("a.u.");
        h_btag_jet2[thr] -> GetXaxis() -> SetTitle("probe jet offline btag (before matching)");
        h_btag_jet2[thr] -> SetStats(0);
        h_btag_jet2[thr] -> SetTitle("");
        h_btag_jet2[thr] -> Scale(1./h_btag_jet2[thr]->GetEntries());
        h_btag_jet2[thr] -> GetYaxis() -> SetRangeUser(0,0.5);
        h_btag_jet2[thr] -> Draw("same hist e");
        h_btag_jet2[thr] -> SetTitle(Form("HLT%s_DoublePFJets%d_PNetBTag_0p11",type.c_str(),thr));
    }
    c_btag_jet2 -> BuildLegend(0.12,0.55,0.60,0.88);
    std::string output_btag_jet2       = input_dir + "/" + input_dir + "_btag_jet2_fh.png";
    if ( atype == 1 ) output_btag_jet2 = input_dir + "/" + input_dir + "_btag_jet2_sl.png";
    c_btag_jet2 -> SaveAs(output_btag_jet2.c_str());

    return 0;
}

int pt_selection(TH1F * h_pt, const float & min_pt, const float & max_pt)
{
    // Get the bin indices corresponding to min_pt and max_pt
    int bin_min = h_pt->GetXaxis()->FindBin(min_pt);
    int bin_max = h_pt->GetXaxis()->FindBin(max_pt);
    for (int bin = 1; bin <= h_pt->GetNbinsX(); ++bin)
    {
        if (bin < bin_min || bin >= bin_max)
        {
            h_pt->SetBinContent(bin, 0);
            h_pt->SetBinError(bin, 0);
        }
    }

    return 0;

}