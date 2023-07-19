#include <typeinfo>
// #include <boost/type_index.hpp>
// #include <boost/core/demangle.hpp>
#include "Rtypes.h"

using rinterface  = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
using RNode = ROOT::RDF::RNode;
using RTH1D = ROOT::RDF::RResultPtr<TH1D>;
using TH1DModel = ROOT::RDF::TH1DModel;

std::map<std::string,int> m_flavour_id_;
std::map<std::string,int> m_flavour_color_;
std::vector<std::string> v_triggers_;
std::map<int,int> m_rebin_;
std::map<int,std::pair<float,float> > m_pt_range_;
std::map<int,std::string> m_pt_range_cut_;
std::map<int,int> m_process_id_;
std::map<int,float> m_probbb_;
std::map<int,std::string> m_probbb_cut_;
std::vector<std::string> v_jettype_;


int init();


int df_analysis()
{
    ROOT::EnableImplicitMT();
    init();
    
    // Define a mutex to protect access to the std::map
    std::mutex mtx;

    bool tag_probbb = true;
    bool jet3_veto = true;
    int workflow = 3;
    bool is_mc = true;

    int color = kRed;
    if ( jet3_veto ) color = kBlue;
    if ( tag_probbb ) color = kBlack;
    std::string s_tag_probbb = "_probbb";
    if ( ! tag_probbb ) s_tag_probbb = "";
    std::string s_jet3_veto = "_j3veto";
    if ( ! jet3_veto ) s_jet3_veto = "";

    std::string era = "*";
    if ( ! is_mc )
    {
        if ( era == "A") color = kMagenta;
        if ( era == "B") color = kCyan;
        if ( era == "C") color = kBlue;
        if ( era == "D") color = kRed;
        if ( era == "*") color = kBlack;
    }

    // open root files
    std::string file_name = Form("nfs/2018/v6/data/wf%d/histograms_2018_v6_data_JetHT_Run2018%s_Jet*_SR_wf%d.root",workflow,era.c_str(),workflow);
    std::string output_name = Form("eff_result%s%s_2018%s",s_jet3_veto.c_str(),s_tag_probbb.c_str(),era.c_str());
    if ( is_mc )
    {
        file_name = Form("nfs/2018/v6/mc/wf%d/histograms_2018_v6_mc_QCD_HT*_Jet*_SR_wf%d.root",workflow,workflow);
        output_name = Form("eff_result%s%s_qcd",s_jet3_veto.c_str(),s_tag_probbb.c_str());
    }
    auto tree_name = "TagAndProbe";
    ROOT::RDataFrame df(tree_name, file_name); 
    // ----

    
    // definitions
    std::function<float(float, float)> pt12_avg = [](float probe_pt, float tag_pt) { return static_cast<float>((probe_pt+tag_pt)/2.f); };
    // ----

    // global filters
    std::function<bool(float, float)> jet3_veto_filter = [&jet3_veto](float jet3_pt, float pt12_avg) { return  !jet3_veto || jet3_pt < 0.3 * pt12_avg; };
    // ----

    // jet pt histograms for each range
    std::map<std::string, std::map<int,RTH1D> > h_pt_jet_range;
    // jet pt histograms for each flavour and range
    std::map<std::string, std::map<int, std::map<std::string,RTH1D> > > h_pt_jet_range_flavour;

    { // mutex guard
        std::lock_guard<std::mutex> lock(mtx);

        auto ndf = df.Define("pt12_avg",pt12_avg,{"probe_pt","tag_pt"})
                    .Filter(jet3_veto_filter, {"jet3_pt", "pt12_avg"})
                    ;
    
        for ( auto & jet_type: v_jettype_) // probe, tag
        {
            for (auto& [range, pt_range_cut]: m_pt_range_cut_)
            {
                int nbins = (int) (m_pt_range_[range].second - m_pt_range_[range].first)/m_rebin_[range];
                std::string histo_cuts = pt_range_cut;
                if ( tag_probbb ) 
                    histo_cuts += Form(" && %s",m_probbb_cut_[range].c_str());
                histo_cuts += Form(" && process_id == %d",m_process_id_[range]);
                if ( jet_type == "probe_match")
                    histo_cuts += " && probe_match == 1";
                
                std::string histo_title = Form("pt_%s_range%d",jet_type.c_str(),range);
                std::string jet_pt_column = Form("%s_pt",jet_type.c_str());
                if ( jet_type == "probe_match" ) jet_pt_column = "probe_pt";
                TH1DModel histo_model = {histo_title.c_str(), histo_title.c_str(), nbins, m_pt_range_[range].first, m_pt_range_[range].second};
                // if ( range > 3)
                    h_pt_jet_range[jet_type].insert({range, ndf.Filter(histo_cuts)
                                                               .Histo1D(histo_model, jet_pt_column,"weight")});
                // else
                //     h_pt_jet_range[jet_type].insert({range, ndf.Filter(histo_cuts)
                //                                            .Filter(jet3_veto_filter, {"jet3_pt", "pt12_avg"})
                //                                            .Histo1D(histo_model, jet_pt_column,"weight")});
                // MC stuff and for tag jets and not for matched probes
                if ( jet_type == "probe_match" || (!is_mc) ) continue;
                for (auto& [flavour, flavour_id]: m_flavour_id_) 
                {
                // if ( range > 3)
                    h_pt_jet_range_flavour[jet_type][range].insert({flavour, ndf.Filter(histo_cuts)
                                                           .Filter(Form("%s_flavour == %d",jet_type.c_str(),flavour_id))
                                                           .Histo1D(histo_model, jet_pt_column,"weight")});
                // else
                //     h_pt_jet_range_flavour[jet_type][range].insert({flavour, ndf.Filter(histo_cuts)
                //                                            .Filter(Form("%s_flavour == %d",jet_type.c_str(),flavour_id))
                //                                            .Filter(jet3_veto_filter, {"jet3_pt", "pt12_avg"})
                //                                            .Histo1D(histo_model, jet_pt_column,"weight")});


                }
            }
        }
    }

    std::map<int,TGraphAsymmErrors*> g_eff_pt_jet_range;
    TMultiGraph * mg_eff_pt_jet = new TMultiGraph();

    std::map<int, std::map<std::string, TGraphAsymmErrors*> > g_frac_pt_probe_range_flavour;
    std::map<int, std::map<std::string, TGraphAsymmErrors*> > g_frac_pt_tag_range_flavour;
    std::map<std::string,TMultiGraph *> mg_frac_pt_probe_flavour;
    std::map<std::string,TMultiGraph *> mg_frac_pt_tag_flavour;
    TMultiGraph * mg_frac_pt_probe = new TMultiGraph();
    TMultiGraph * mg_frac_pt_tag   = new TMultiGraph();

    TLegend * leg_frac_probe = new TLegend(0.75,0.65,0.85,0.85,NULL,"brNDC");
    TLegend * leg_frac_tag   = new TLegend(0.15,0.35,0.25,0.55,NULL,"brNDC");



    for (auto& [flavour, flavour_color]: m_flavour_color_)
    {
        mg_frac_pt_probe_flavour[flavour] = new TMultiGraph();
        mg_frac_pt_tag_flavour[flavour]   = new TMultiGraph();
    }

    {   // mutex guard
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& [range, pt_range]: m_pt_range_)
        {
            g_eff_pt_jet_range[range] = new TGraphAsymmErrors(h_pt_jet_range["probe_match"][range].GetPtr(),h_pt_jet_range["probe"][range].GetPtr(),"cl=0.683 b(1,1) mode");
            g_eff_pt_jet_range[range]->SetMarkerStyle(20);
            g_eff_pt_jet_range[range]->SetMarkerColor(color);
            g_eff_pt_jet_range[range]->SetLineColor(color);
            
            mg_eff_pt_jet->Add(g_eff_pt_jet_range[range],"p");

            if ( ! is_mc ) continue;
            for (auto& [flavour, flavour_color]: m_flavour_color_)
            {
                g_frac_pt_probe_range_flavour[range][flavour] = new TGraphAsymmErrors(h_pt_jet_range_flavour["probe"][range][flavour].GetPtr(),h_pt_jet_range["probe"][range].GetPtr(),"cl=0.683 b(1,1) mode");
                g_frac_pt_probe_range_flavour[range][flavour] -> SetMarkerStyle(20);
                g_frac_pt_probe_range_flavour[range][flavour] -> SetMarkerColor(flavour_color);
                g_frac_pt_probe_range_flavour[range][flavour] -> SetLineColor(flavour_color);
                if ( range == 1 )
                {
                    leg_frac_probe -> AddEntry(g_frac_pt_probe_range_flavour[range][flavour],flavour.c_str(),"p");
                }
                mg_frac_pt_probe_flavour[flavour]->Add(g_frac_pt_probe_range_flavour[range][flavour],"P");
                mg_frac_pt_probe->Add(g_frac_pt_probe_range_flavour[range][flavour],"P");
                g_frac_pt_tag_range_flavour[range][flavour]   = new TGraphAsymmErrors(h_pt_jet_range_flavour["tag"][range][flavour].GetPtr(),h_pt_jet_range["tag"][range].GetPtr(),"cl=0.683 b(1,1) mode");
                g_frac_pt_tag_range_flavour[range][flavour] -> SetMarkerStyle(20);
                g_frac_pt_tag_range_flavour[range][flavour] -> SetMarkerColor(flavour_color);
                g_frac_pt_tag_range_flavour[range][flavour] -> SetLineColor(flavour_color);
                if ( range == 1 )
                {
                    leg_frac_tag -> AddEntry(g_frac_pt_tag_range_flavour[range][flavour],flavour.c_str(),"p");
                }

                mg_frac_pt_tag_flavour[flavour]->Add(g_frac_pt_tag_range_flavour[range][flavour],"P");
                mg_frac_pt_tag->Add(g_frac_pt_tag_range_flavour[range][flavour],"P");

            }

            
        }

    }
    
    {
        std::lock_guard<std::mutex> lock(mtx);
        TCanvas * c_eff = new TCanvas("c_eff","efficiency");
        mg_eff_pt_jet->SetName("online_btag_efficiency");
        mg_eff_pt_jet->GetHistogram()->GetXaxis()->SetTitle("jet p_{T} [GeV]");
        mg_eff_pt_jet->GetHistogram()->GetYaxis()->SetTitle("rel. online btag efficiency");
        mg_eff_pt_jet->Draw("a");
        mg_eff_pt_jet -> GetHistogram() -> GetYaxis() -> SetRangeUser(0.,1.);
        c_eff->Modified();
        c_eff->Update();
        c_eff->SaveAs(Form("df_analysis_%s.png",output_name.c_str()));
        TFile * f_out = new TFile(Form("df_analysis_%s.root",output_name.c_str()),"recreate");
        mg_eff_pt_jet->Write();

        if ( is_mc )
        {
            TCanvas * c_frac_probe = new TCanvas("c_frac_probe","");
            mg_frac_pt_probe->SetName("mg_frac_pt_probe");
            mg_frac_pt_probe-> GetHistogram() -> GetYaxis() -> SetRangeUser(0.,1.);
            mg_frac_pt_probe-> GetHistogram() -> GetYaxis() -> SetTitle("flavour fraction of probe jets");
            mg_frac_pt_probe-> GetHistogram() -> GetXaxis() -> SetTitle("probe jet p_{T} [GeV]");
            mg_frac_pt_probe->Draw("AP");
            leg_frac_probe->Draw();
            mg_frac_pt_probe->Write();
            c_frac_probe -> SaveAs(Form("df_analysis_frac_probe_%s.png",output_name.c_str()));
            TCanvas * c_frac_tag = new TCanvas("c_frac_tag","");
            mg_frac_pt_tag->SetName("mg_frac_pt_tag");
            mg_frac_pt_tag-> GetHistogram() -> GetYaxis() -> SetRangeUser(0.,1.);
            mg_frac_pt_tag-> GetHistogram() -> GetYaxis() -> SetTitle("flavour fraction of tag jets");
            mg_frac_pt_tag-> GetHistogram() -> GetXaxis() -> SetTitle("tag jet p_{T} [GeV]");
            mg_frac_pt_tag->Draw("AP");
            leg_frac_tag->Draw();
            mg_frac_pt_tag->Write();
            c_frac_tag -> SaveAs(Form("df_analysis_frac_tag_%s.png",output_name.c_str()));


        }
        f_out->Close();

    }
        

 
    return 0;
}


int init()
{
    // colors
    m_flavour_color_["b"]    = kRed;
    m_flavour_color_["bb"]   = kMagenta;
    m_flavour_color_["c"]    = kBlue;
    m_flavour_color_["cc"]   = kCyan;
    m_flavour_color_["udsg"] = kBlack;
    m_flavour_id_["b"]    = 5;
    m_flavour_id_["bb"]   = 55;
    m_flavour_id_["c"]    = 4;
    m_flavour_id_["cc"]   = 44;
    m_flavour_id_["udsg"] = 0;
    // triggers aliases
    v_triggers_.push_back("Jet40" );
    v_triggers_.push_back("Jet100");
    v_triggers_.push_back("Jet200");
    v_triggers_.push_back("Jet350");
    // m_pt_range_[1] = std::make_pair(30.,40.);
    m_pt_range_[1] = std::make_pair(40.,80.);
    m_pt_range_[2] = std::make_pair(80.,130.);
    m_pt_range_[3] = std::make_pair(130.,250.);
    m_pt_range_[4] = std::make_pair(250.,550.);
    m_pt_range_[5] = std::make_pair(550.,800.);
    // m_pt_range_[6] = std::make_pair(800.,1000.);

    // binning for each range
    // m_rebin_[1] = 10;
    // m_rebin_[2] = 10;
    // m_rebin_[3] = 10;
    // m_rebin_[4] = 20;
    // m_rebin_[5] = 50;
    // m_rebin_[6] = 100;
    // m_rebin_[7] = 200;
    // m_rebin_[1] = 10;
    m_rebin_[1] = 20;
    m_rebin_[2] = 25;
    m_rebin_[3] = 60;
    m_rebin_[4] = 100;
    m_rebin_[5] = 250;
    m_rebin_[6] = 200;
        
    // m_probbb_[1] = 0.02;
    m_probbb_[1] = 0.02;
    m_probbb_[2] = 0.02;
    m_probbb_[3] = 0.002;
    m_probbb_[4] = 0.002;
    m_probbb_[5] = 0.002;
    m_probbb_[6] = 0.002;

    // m_process_id_[1] = 1;
    m_process_id_[1] = 1;
    m_process_id_[2] = 1;
    m_process_id_[3] = 2;
    m_process_id_[4] = 3;
    m_process_id_[5] = 3;
    m_process_id_[6] = 4;

    for (auto& [bin, range]: m_pt_range_)
    {
        m_pt_range_cut_[bin] = Form("probe_pt>=%.2f  && probe_pt<%.2f", range.first, range.second);
        m_probbb_cut_[bin] = Form("tag_prob_bb<%.5f",m_probbb_[bin]); // tag prob bb
    }


    
    // jet type
    v_jettype_.push_back("tag");
    v_jettype_.push_back("probe");
    v_jettype_.push_back("probe_match");
    
    
    return 0;
}
