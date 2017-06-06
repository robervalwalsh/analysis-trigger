import os

triggers_sl = []

triggers_sl.append('HLT_1CaloJets30_Muon12_1PFJets40_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_2CaloBTagCSV084_2PFJets40_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets40_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets50_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets60_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_2CaloBTagCSV092_2PFJets70_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets40_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets100_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets200_psw')
triggers_sl.append('HLT_2CaloJets30_Muon12_1CaloBTagCSV092_2PFJets350_psw')


triggers_ah = []
triggers_ah.append('HLT_2CaloJets100_2CaloBTagCSV084_2PFJets100_psw')
triggers_ah.append('HLT_2CaloJets100_2CaloBTagCSV092_2PFJets100_psw')
triggers_ah.append('HLT_2CaloJets100_2CaloBTagCSV092_2PFJets110_psw')
triggers_ah.append('HLT_2CaloJets100_2CaloBTagCSV092_2PFJets120_psw')
triggers_ah.append('HLT_2CaloJets100_2CaloBTagCSV092_2PFJets130_psw')
triggers_ah.append('HLT_2CaloJets100_1CaloBTagCSV092_2PFJets100_psw')
triggers_ah.append('HLT_2CaloJets100_1CaloBTagCSV092_2PFJets140_psw')
triggers_ah.append('HLT_2CaloJets100_1CaloBTagCSV092_2PFJets200_psw')
triggers_ah.append('HLT_2CaloJets100_1CaloBTagCSV092_2PFJets350_psw')

triggers_jh = []
triggers_jh.append('HLT_2CaloJets30_1CaloBTagCSV092_2PFJets40_psw')


for trigger in triggers_sl:
   command = 'root -l -q -b HLTRatesMC.C\(\\"'+trigger+'\\"\)'
   print command
   os.system(command)
   
for i,trigger in enumerate(triggers_ah):
   command = 'root -l -q -b HLTRatesMC.C\(\\"'+trigger+'\\"\)'
   print command
   os.system(command)
   
for trigger in triggers_jh:
   command = 'root -l -q -b HLTRatesMC.C\(\\"'+trigger+'\\"\)'
   print command
   os.system(command)
