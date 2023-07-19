
TGraphAsymmErrors* ratioGraph(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2);
TMultiGraph* ratioMultiGraph(TMultiGraph* mg1, TMultiGraph* mg2);
TGraphAsymmErrors* combineGraphPoints(TMultiGraph* mg);

int eff_plots()
{   
    // std::vector<std::string> v_era = {"A","B","C","D","*"};
    std::vector<std::string> v_era = {"qcd","*"};
    std::vector<int> v_color = {kRed,kBlue,kMagenta,kCyan,kBlack};
    std::vector<TFile *> v_file;
    std::vector<TMultiGraph *> v_mg_eff;

    TMultiGraph * mg_eff = new TMultiGraph();

    bool j3veto = true;
    bool probbb = true;

    std::string name_sel = "";
    if ( j3veto )
    {
        name_sel = "j3veto";
        if ( probbb ) name_sel+= "_probbb";
    }

    for ( size_t i = 0; i < v_era.size() ; ++i )
    {
        auto era = v_era[i];
        auto color = v_color[i];
        std::string filename        = Form("df_analysis_eff_result_%s_2018%s.root",name_sel.c_str(),era.c_str());
        if ( era == "qcd") filename = Form("df_analysis_eff_result_%s_qcd.root",name_sel.c_str());
        v_file.push_back(new TFile(filename.c_str(),"old"));
        v_mg_eff.push_back((TMultiGraph*) v_file.back()->Get("online_btag_efficiency"));

        TList* graphList = v_mg_eff.back()->GetListOfGraphs();
        if (graphList) {
            Int_t numGraphs = graphList->GetSize();
            for (Int_t i = 0; i < numGraphs; ++i) {
                TGraphAsymmErrors* graph = dynamic_cast<TGraphAsymmErrors*>(graphList->At(i));
                if (graph) {
                    graph -> SetMarkerColor(color);
                    graph -> SetLineColor(color);
                }
            }
        }        
        mg_eff->Add(v_mg_eff.back(),"p");
    }

    TCanvas * c1 = new TCanvas("c1","");
    
    mg_eff->GetHistogram()->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    mg_eff->GetHistogram()->GetYaxis()->SetTitle("rel. online btag efficiency");
    mg_eff->GetHistogram()->GetYaxis()->SetRangeUser(0,1);

    mg_eff->Draw("a");

    if ( v_era[0] != "qcd")
        c1 -> SaveAs(Form("eff_data_compare_%s_2018.png",name_sel.c_str()));
    else
        c1 -> SaveAs(Form("eff_data_mc_%s_2018.png",name_sel.c_str()));

    TMultiGraph * mg_sf = ratioMultiGraph(v_mg_eff[1],v_mg_eff[0]);

    TGraphAsymmErrors* g_sf = combineGraphPoints(mg_sf);
    g_sf -> SetMarkerStyle(20);
    g_sf -> SetName("scale_factors_graph");

    TF1 *fpol = new TF1("fpol", "pol1", 30, 1000);
    fpol -> SetLineColor(kBlack);
    fpol -> SetLineStyle(2);
    auto fit_result = g_sf->Fit(fpol, "QS");
    fpol -> FixParameter(0,fpol->GetParameter(0));
    fit_result = g_sf->Fit(fpol, "QS");

    //Create a histo with errors to store the intervals
    int ngr = g_sf->GetN();
    TH1D *hint = new TH1D("hint","", 760, 40, 800);
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint,0.68);

    int n_fit_sf = hint->GetNbinsX();
    TGraphErrors * g_fit_sf = new TGraphErrors(n_fit_sf);
    for ( int i = 0; i < hint->GetNbinsX(); i++ )
    {
        g_fit_sf -> SetPoint(i,hint->GetBinCenter(i+1),hint->GetBinContent(i+1));
        g_fit_sf -> SetPointError(i,0.,hint->GetBinError(i+1));
        // std::cout << "bin " << i+1 << ":  " << hint->GetBinCenter(i+1) << "  " << hint->GetBinContent(i+1) << " ± " << hint->GetBinError(i+1) << std::endl;
    }

    hint->SetFillColor(2);
    hint->SetFillStyle(3001);
    hint->GetYaxis()->SetRangeUser(0.3,1.7);
    hint->GetYaxis()->SetTitle("rel. online btag scale factor");
    hint->GetXaxis()->SetTitle("jet p_{T} [GeV]");
    hint->SetStats(false);
    hint->Draw("e3");
    g_sf->Draw("psame");

    c1 ->SaveAs(Form("online_btag_SF_%s_2018.png",name_sel.c_str()));

    std::cout << fit_result -> Chi2()/fit_result->Ndf() << std::endl;

    fit_result ->Print("V");

    g_fit_sf -> SetName("fit_result_graph_error_1sigma");
    TFile * f_fit_out = new TFile(Form("online_btag_SF_%s_2018.root",name_sel.c_str()),"recreate");
    fit_result -> Write();
    fpol->Write();
    g_fit_sf -> Write();
    g_sf->Write();
    f_fit_out -> Close();
    return 0;
}


TGraphAsymmErrors* ratioGraph(TGraphAsymmErrors* graph1, TGraphAsymmErrors* graph2) {
    int numPoints1 = graph1->GetN();
    int numPoints2 = graph2->GetN();

    if (numPoints1 != numPoints2) {
        // Error handling: graphs should have the same number of points
        return nullptr;
    }

    TGraphAsymmErrors* ratio_graph = new TGraphAsymmErrors(numPoints1);

    for (int i = 0; i < numPoints1; ++i) {
        // Retrieve the x and y values of the first graph
        double x1, y1;
        graph1->GetPoint(i, x1, y1);

        // Retrieve the x and y values of the second graph
        double x2, y2;
        graph2->GetPoint(i, x2, y2);

        // Calculate the ratio of y values
        double ratio = y1 / y2;

        // Calculate the error on the ratio using error propagation assuming uncorrelated uncertainties
        double err1Up = graph1->GetErrorYhigh(i);
        double err1Down = graph1->GetErrorYlow(i);
        double err2Up = graph2->GetErrorYhigh(i);
        double err2Down = graph2->GetErrorYlow(i);

        double errRatioUp = std::sqrt(std::pow(err1Up / y2, 2) + std::pow(y1 * err2Up / (y2 * y2), 2));
        double errRatioDown = std::sqrt(std::pow(err1Down / y2, 2) + std::pow(y1 * err2Down / (y2 * y2), 2));

        // Set the point and error values in the ratio graph
        ratio_graph->SetPoint(i, x1, ratio);
        ratio_graph->SetPointError(i, 0., 0., errRatioDown, errRatioUp);
    }

    return ratio_graph;
}

TMultiGraph* ratioMultiGraph(TMultiGraph* mg1, TMultiGraph* mg2) {
    if (mg1->GetListOfGraphs()->GetSize() != mg2->GetListOfGraphs()->GetSize()) {
        // Error handling: multi-graphs should have the same number of graphs
        return nullptr;
    }

    TMultiGraph* ratioMg = new TMultiGraph();

    TIter nextGraph1(mg1->GetListOfGraphs());
    TIter nextGraph2(mg2->GetListOfGraphs());

    TGraphAsymmErrors* graph1;
    TGraphAsymmErrors* graph2;

    while ((graph1 = dynamic_cast<TGraphAsymmErrors*>(nextGraph1()))) {
        graph2 = dynamic_cast<TGraphAsymmErrors*>(nextGraph2());

        TGraphAsymmErrors* ratio_graph = ratioGraph(graph1, graph2);

        if (ratio_graph) {
            ratioMg->Add(ratio_graph);
        } else {
            // Error handling: ratio calculation failed
            delete ratioMg;
            return nullptr;
        }
    }

    return ratioMg;
}

TGraphAsymmErrors* combineGraphPoints(TMultiGraph* mg) {
    TGraphAsymmErrors* combinedGraph = new TGraphAsymmErrors();
    
    TIter nextGraph(mg->GetListOfGraphs());
    TGraphAsymmErrors* graph;

    while ((graph = dynamic_cast<TGraphAsymmErrors*>(nextGraph()))) {
        Int_t nPoints = graph->GetN();
        
        Double_t* x = graph->GetX();
        Double_t* y = graph->GetY();
        Double_t* exl = graph->GetEXlow();
        Double_t* exh = graph->GetEXhigh();
        Double_t* eyl = graph->GetEYlow();
        Double_t* eyh = graph->GetEYhigh();
        
        for (Int_t i = 0; i < nPoints; ++i) {
            Double_t xValue = x[i];
            Double_t yValue = y[i];
            Double_t exlValue = exl[i];
            Double_t exhValue = exh[i];
            Double_t eylValue = eyl[i];
            Double_t eyhValue = eyh[i];
            
            combinedGraph->SetPoint(combinedGraph->GetN(), xValue, yValue);
            combinedGraph->SetPointError(combinedGraph->GetN() - 1, exlValue, exhValue, eylValue, eyhValue);
        }
    }
    
    return combinedGraph;
}
