#include <typeinfo>
#include <boost/type_index.hpp>
#include <boost/core/demangle.hpp>
#include "Rtypes.h"

using rinterface  = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;
using rnode = ROOT::RDF::RNode;

std::map<std::string,int> m_flavour_id_;
std::map<std::string,int> m_flavour_color_;
std::vector<std::string> v_triggers_;
std::map<int,int> m_rebin_;
std::map<int,std::pair<float,float> > m_pt_range_;
std::map<int,std::string> m_pt_range_cut_;
std::map<int,float> m_probbb_;
std::map<int,std::string> m_probbb_cut_;
std::vector<std::string> v_jettype_;


int init();


int tree_analysis()
{
    init();

    bool tag_probbb = true;
    int workflow = 2;

    auto file_name = Form("nfs/2018/v6/mc/save/wf%d/histograms_2018_v6_mc_QCD_HT*_Jet*_SR_wf%d.root",workflow,workflow);
    auto tree_name = "TagAndProbe";
    ROOT::RDataFrame df(tree_name, file_name); 

    //~ std::map<int,rnode> m_rdf;
    //~ std::map<int,TH1D*> h_probe_pt_range_all;
    //~ std::map<int,std::map<std::string,TH1D*> > h_probe_pt_range_flavour;
    //~ std::map<int,std::map<std::string,TGraphAsymmErrors*> > g_probe_range_flavour_fraction;
    //~ TMultiGraph *mg_probe_flavour_fraction = new TMultiGraph();
    
    std::map<int,rnode> m_rdf;
    std::map<std::string, std::map<int,TH1D*> > h_pt_jet_range_all;
    
    // flavour dependency
    std::map<std::string, std::map<int,std::map<std::string,TH1D*> > > h_pt_jet_range_flavour;
    std::map<std::string, std::map<int,std::map<std::string,TGraphAsymmErrors*> > > g_fraction_jet_range_flavour;
    std::map<std::string, TMultiGraph*> mg_fraction_jet_flavour;
    mg_fraction_jet_flavour["probe"] = new TMultiGraph();
    mg_fraction_jet_flavour["tag"]   = new TMultiGraph();
    
    
    

    for (auto& [range, pt_cut]: m_pt_range_cut_)
    {
        // one rnode per pt range
        std::string cut = pt_cut;
        if ( tag_probbb ) 
            cut += Form(" && %s",m_probbb_cut_[range].c_str());
        m_rdf.insert({range,rnode(df.Filter(cut))});
        auto bin_size = (int) (m_pt_range_[range].second - m_pt_range_[range].first)/m_rebin_[range];
        for ( auto & jet_type: v_jettype_)
        {
            std::string histo_title = Form("pt_%s_range%d",jet_type.c_str(),range);
            h_pt_jet_range_all[jet_type][range] = new TH1D(m_rdf.at(range).Histo1D({histo_title.c_str(), histo_title.c_str(), bin_size, m_pt_range_[range].first, m_pt_range_[range].second}, Form("%s_pt",jet_type.c_str()),"weight").GetValue());
        
            for ( auto & [flavour, id]: m_flavour_id_)
            {
                std::string flavour_filter = Form("%s_flavour==%d",jet_type.c_str(),id);
                std::string histo_title_flavour = Form("pt_%s_range%d_%s",jet_type.c_str(),range,flavour.c_str());
                h_pt_jet_range_flavour[jet_type][range][flavour] = new TH1D(m_rdf.at(range).Filter(flavour_filter).Histo1D({histo_title_flavour.c_str(), histo_title_flavour.c_str(), bin_size, m_pt_range_[range].first, m_pt_range_[range].second}, Form("%s_pt",jet_type.c_str()),"weight").GetValue());
            
                g_fraction_jet_range_flavour[jet_type][range][flavour] = new TGraphAsymmErrors(h_pt_jet_range_flavour[jet_type][range][flavour],h_pt_jet_range_all[jet_type][range],"cl=0.683 b(1,1) mode");
                g_fraction_jet_range_flavour[jet_type][range][flavour]->SetMarkerColor(m_flavour_color_[flavour]);
                g_fraction_jet_range_flavour[jet_type][range][flavour]->SetMarkerStyle(20);
                g_fraction_jet_range_flavour[jet_type][range][flavour]->SetLineColor(m_flavour_color_[flavour]);
                
                mg_fraction_jet_flavour[jet_type]->Add(g_fraction_jet_range_flavour[jet_type][range][flavour],"p");
            
                //~ if ( flavour == "b" )
                //~ {
                    //~ auto n = g_fraction_jet_range_flavour[jet_type][range][flavour] -> GetN();
                    //~ auto x = g_fraction_jet_range_flavour[jet_type][range][flavour]-> GetX();
                    //~ auto y = g_fraction_jet_range_flavour[jet_type][range][flavour] -> GetY();
                    //~ for ( int i = 0; i < n; ++i)
                        //~ cout << x[i] << "   " << y[i] << endl;
                //~ }
                
            }
        }
    }
    
    std::map<std::string, TCanvas*> c_fraction;
    for (auto & jet_type : v_jettype_)
    {
        std::string c_fraction_title = Form("c_fraction_%s",jet_type.c_str());
        c_fraction[jet_type] = new TCanvas(c_fraction_title.c_str(),c_fraction_title.c_str());
        mg_fraction_jet_flavour[jet_type] -> Draw("a");
        mg_fraction_jet_flavour[jet_type] -> GetHistogram() -> GetYaxis() -> SetRangeUser(0.,1.);
        c_fraction[jet_type]->Modified();
        c_fraction[jet_type]->Update();
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
    m_pt_range_[1] = std::make_pair(30.,40.);
    m_pt_range_[2] = std::make_pair(40.,80.);
    m_pt_range_[3] = std::make_pair(80.,130.);
    m_pt_range_[4] = std::make_pair(130.,250.);
    m_pt_range_[5] = std::make_pair(250.,600.);
    m_pt_range_[6] = std::make_pair(600.,800.);
    m_pt_range_[7] = std::make_pair(800.,1000.);
    m_probbb_[1] = 0.005;
    m_probbb_[2] = 0.02;
    m_probbb_[3] = 0.02;
    m_probbb_[4] = 0.002;
    m_probbb_[5] = 0.002;
    m_probbb_[6] = 0.002;
    m_probbb_[7] = 0.002;
    for (auto& [bin, range]: m_pt_range_)
    {
        m_pt_range_cut_[bin] = Form("probe_pt>=%.2f  && probe_pt<%.2f", range.first, range.second);
        m_probbb_cut_[bin] = Form("tag_prob_bb<%.5f",m_probbb_[bin]); // tag prob bb
    }

    // binning for each range
    m_rebin_[1] = 10;
    m_rebin_[2] = 10;
    m_rebin_[3] = 10;
    m_rebin_[4] = 20;
    m_rebin_[5] = 50;
    m_rebin_[6] = 100;
    m_rebin_[7] = 200;
    
    // jet type
    v_jettype_.push_back("probe");
    v_jettype_.push_back("tag");
    
    
    return 0;
}
