#include "Analysis/Tools/interface/Analyser.h"

using namespace analysis;
using namespace analysis::tools;

int main(int argc, char ** argv)
{
    TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

    // construct your analyser
    Analyser analyser(argc,argv);
    
    analyser.jetHistograms("tandp_selection");
    analyser.jetHistograms("tandp_probe_matched");
    

    for ( int i = 0 ; i < analyser.nEvents() ; ++i )
    {
        if ( ! analyser.event(i)              )  continue;    // read event, run selection/json
        if ( ! analyser.triggerSelection ()   )  continue;    // trigger
        if ( ! analyser.preselection()        )  continue;    // preselection - particle id's
        if ( ! analyser.jetCorrections()      )  continue;    // jet corrections
        
        // select events according to # jets criteria in the config file
        if ( ! analyser.selectionNJets()      )  continue;
        // select events according to criteria for the leading and pr2nd leading jets kinematics 
        if ( ! analyser.selectionJet(1)       )  continue;
        if ( ! analyser.selectionJet(2)       )  continue;
        // select events according to delta_phi criteria in the config file
        if ( ! analyser.selectionJetDphi(1,2) )  continue;
        
        ///////////////////////////////////
        // TO DO: veto events with a third jet
        ///////////////////////////////////
        
        // Assuming the tag to be the leading jet,
        // select events according to matching the leading jet to online objects defined in the config file
        if ( ! analyser.onlineJetMatching(1) )  continue;
        
        // Fill predefined histograms for the event selection
        // In the directory 'selection', the histograms for jet 2 correspond probe jets
        analyser.fillJetHistograms("tandp_selection");
        
        // select events where the 2nd jet has a match to the online objects defined in the configuration
        if ( ! analyser.onlineJetMatching(2) )  continue;
        // Fill predefined histograms for the event selection after the probe jets are matched to online objects
        // In the directory 'selection_probe_matched', the histograms for jet 2 are those for the probe jets after matching to online
        analyser.fillJetHistograms("tandp_probe_matched");
        
        
        
        // 
              
    }
} // end main
     



