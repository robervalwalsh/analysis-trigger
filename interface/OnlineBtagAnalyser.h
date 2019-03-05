#ifndef Analysis_Trigger_OnlineBtagAnalyser_h
#define Analysis_Trigger_OnlineBtagAnalyser_h 1

// -*- C++ -*-
//
// Package:    Analysis/Trigger
// Class:      Analysis
// 
/**\class Analysis OnlineBtagAnalyser.cc Analysis/Trigger/src/OnlineBtagAnalyser.cc

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
#include <memory>
#include <vector>
#include <string>
// 
// user include files

#include "Analysis/Tools/interface/Analyser.h"

//
// class declaration
//

namespace analysis {
   namespace trigger {

      class OnlineBtagAnalyser : public analysis::tools::Analyser {
         
         public:
            OnlineBtagAnalyser();
            OnlineBtagAnalyser(int argc, char ** argv);
           ~OnlineBtagAnalyser();
           
            // ----------member data ---------------------------
         protected:
               
         private:
               
         public:
           virtual bool muonJet(const bool & swap = true);

      };
   }
}

#endif  // Analysis_Trigger_OnlineBtagAnalyser_h
