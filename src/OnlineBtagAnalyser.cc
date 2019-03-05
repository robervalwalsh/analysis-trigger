/**\class Trigger OnlineBtagAnalyser.cc Analysis/Trigger/src/OnlineBtagAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh
//         Created:  Mon, 05 Mar 2019 18:30:00 GMT
//
//

// system include files
#include <iostream>
// 
// user include files
#include "Analysis/Trigger/interface/OnlineBtagAnalyser.h"

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::trigger;

//
// constructors and destructor
//
OnlineBtagAnalyser::OnlineBtagAnalyser()
{
}

OnlineBtagAnalyser::OnlineBtagAnalyser(int argc, char ** argv) : BaseAnalyser(argc,argv), Analyser(argc,argv)
{
}

OnlineBtagAnalyser::~OnlineBtagAnalyser()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
}


//
// member functions
//
// ------------ method called for each event  ------------


bool OnlineBtagAnalyser::muonJet(const bool & swap)
{
   // jet with ranking 1 is the muon jet, swap with jet 2 in case it has the leading muon 
   int r1 = 1;
   int r2 = 2;
   int j1 = r1-1;
   int j2 = r2-1;
   ++ cutflow_;
   if ( std::string(h1_["cutflow"] -> GetXaxis()-> GetBinLabel(cutflow_+1)) == "" ) 
   {
      if ( swap ) h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Online Btag T&P Semileptonic: Jet-muon association -> Muon-Jet index 1 (probe jet)");
      else        h1_["cutflow"] -> GetXaxis()-> SetBinLabel(cutflow_+1,"Online Btag T&P Semileptonic: Jet-muon association");
   }
   
   auto jet1 = selectedJets_[j1];
   jet1 -> addMuon(selectedMuons_);
   auto jet2 = selectedJets_[j2];
   jet2 -> addMuon(selectedMuons_);
   
   if ( ! (jet1 -> muon() || jet2 -> muon()) ) return false;
   if ( !  jet1 -> muon() && swap ) this->jetSwap(r1,r2);
   
   h1_["cutflow"] -> Fill(cutflow_,weight_);
   return true;
   
}

