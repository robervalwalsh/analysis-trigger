int pt_plots()
{   
    int workflow = 1;
    TFile * f[4];
    TH1F * h_pt[4];
    std::string trigger[4] = {"Jet40", "Jet100", "Jet200", "Jet350"};
    int color[4] = {kRed, kBlue, kMagenta, kBlack};
    std::string dir = Form("/nfs/dust/cms/user/walsh/issues/cms/11/Analysis/Trigger/test/run2_ul/2018/v6/data/wf%d/histograms/no_pt_cut",workflow);
    dir = Form("/Users/walsh/nfs/dust/cms/user/walsh/issues/cms/11/Analysis/Trigger/test/run2_ul/2018/v6/data/wf%d/histograms/no_pt_cut",workflow);

    TCanvas * c1 = new TCanvas("c1","",1000,1000);
    c1->SetLogy();
    for ( int i = 0; i < 4; ++i)
    {
        f[i] = new TFile(Form("%s/histograms_2018_v6_data_JetHT_Run2018D_%s_SR_wf%d.root",dir.c_str(),trigger[i].c_str(),workflow),"old");
        h_pt[i] = (TH1F*) f[i] -> Get("probe/pt_jet1");
        h_pt[i] -> SetTitle("");
        h_pt[i] -> SetLineColor(color[i]);
        h_pt[i] -> SetLineWidth(2);
        h_pt[i] -> SetMarkerStyle(20);
        h_pt[i] -> SetMarkerColor(color[i]);
        h_pt[i] -> Rebin(10);
        h_pt[i] -> GetXaxis()->SetRangeUser(0,800);
        h_pt[i] -> Draw("same");
    }
    
    c1 ->SaveAs(Form("tag_and_probe_wf%d_pt_no_cut.png",workflow));

    // TH1F * h_pt_sum = (TH1F*) h_pt[0]->Clone("h_pt_sum");
    // h_pt_sum->Add(h_pt[1]);
    // h_pt_sum->Add(h_pt[2]);
    // h_pt_sum->Add(h_pt[3]);

    // TCanvas * c2 = new TCanvas("c2","",1000,1000);
    // c2->SetLogy();
    // h_pt_sum->Draw();

    return 0;
}