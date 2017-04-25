#!/bin/csh -f

# DATA
# set file = "../test/80x/BTagMu20/mssmhbb_triggers_data_BTagMu20.root"
# set trigger =    "HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger =    "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger =    "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger =    "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# 
# set file = "../test/80x/BTagCSV/mssmhbb_triggers_data_BTagCSV.root"
# set trigger = "HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)
# set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5_psw"
# root -l -q -b HLTRatesData.C\(\"$trigger\",\"$file\"\)


# MC
set trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50e40_Eta2p1_dEta1p5_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)

set trigger = "HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
set trigger = "HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_Eta2p1_dEta1p5_psw"
root -l -q -b HLTRatesMC.C\(\"$trigger\"\)
