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



files = []
files.append('../test/80x/BTagMu20/mssmhbb_triggers_data_BTagMu20_control_backup_v5.root')
files.append('../test/80x/BTagCSV/mssmhbb_triggers_data_BTagCSV_control_backup_v4.root')
files.append('../test/80x/BTagCSVControl/mssmhbb_triggers_data_BTagCSVControl_control_backup_v4.root')
files.append('../test/80x/JetHT/mssmhbb_triggers_data_JetHT_control_backup_v5.root')

for trigger in triggers_sl:
   command = 'root -l -q -b HLTRatesData.C\(\\"'+trigger+'\\",\\"'+files[0]+'\\"\)'
   print command
   os.system(command)
   
for i,trigger in enumerate(triggers_ah):
   if i <= 4:
      command = 'root -l -q -b HLTRatesData.C\(\\"'+trigger+'\\",\\"'+files[1]+'\\"\)'
   else:
      command = 'root -l -q -b HLTRatesData.C\(\\"'+trigger+'\\",\\"'+files[2]+'\\"\)'
   
   print command
   os.system(command)
   
for trigger in triggers_jh:
   command = 'root -l -q -b HLTRatesData.C\(\\"'+trigger+'\\",\\"'+files[3]+'\\"\)'
   print command
   os.system(command)
   
