#!/usr/bin/env python3

import os, glob, sys
import pandas as pd

# todo: read csv and template cfg from the command line
# todo: as well as data or qcd

current_dir = os.getcwd()
nfs_dir = os.path.realpath('nfs')

def create_config(alias, parameters, template):
    config = template.replace('template.cfg', f'{alias}.cfg')
    # Read in the file
    with open(template, 'r') as file :
        filedata = file.read()
    # Replace the target string
    for par, value in parameters.items():
        filedata = filedata.replace(par, str(value))
    # Write the file out again
    with open(config, 'w') as file:
        file.write(filedata)
    return config

def get_parameters(df,index):
    headers = list(df)
    alias = df.iloc[index,0]
    parameters = {}
    for idx,h in enumerate(headers):
        if idx == 0: continue
        parameters[h] = df.iloc[index, idx]
    return alias,parameters

def get_files_list(year,version,type,benriched = False):
    cmssw_src = os.path.join(os.getenv('CMSSW_BASE'),'src')
    data_dir = os.path.join(cmssw_src,f'Analysis/Tools/data/ntuples/{year}/{version}/{type}')
    files_list = []
    if type == 'mc':
        files_list_bEnriched = glob.glob(f'{data_dir}/QCD_bEnriched_HT*_rootFileList.txt')
        files_list_BGenFilter = glob.glob(f'{data_dir}/QCD_HT*_BGenFilter_rootFileList.txt')
        files_list_QCD = glob.glob(f'{data_dir}/QCD_HT*_rootFileList.txt')
        files_list_QCD = [f for f in files_list_QCD if 'BGenFilter' not in f]
        files_list = files_list_QCD
        if benriched:
            files_list = files_list_bEnriched+files_list_BGenFilter
    if type == 'data':
        files_list = glob.glob(f'{data_dir}/JetHT*_rootFileList.txt')
    return files_list

configs = pd.read_csv('configs_2018.csv')
template = 'tagandprobe_fh_2018ul_template.cfg'
year = '2018'
version = 'v6'
type = 'data'
workflow = 4
workflow_title = 'Tag&Probe online btag'
if workflow == 2:
    workflow_title += ' w/ 3rd jet veto'

cmssw_src = os.path.join(os.getenv('CMSSW_BASE'),'src')
json = f'{cmssw_src}/Analysis/Tools/data/calibrations/{year}/certified/Cert_Run2018ABCD-DeepCSV_13TeV_Legacy2018_Collisions18_JSON.txt'

work_dir = f'{nfs_dir}/{year}/{version}/{type}/wf{workflow}'
os.makedirs(work_dir, exist_ok=True)

config_files = {}
for i in range(len(configs)):
    alias,parameters = get_parameters(configs,i)
    config_file = create_config(alias,parameters,template)
    config_files[alias] = config_file
    try:
        os.symlink(f'{current_dir}/{config_file}',f'{work_dir}/{config_file}')
    except FileExistsError:
        pass

# print(*files_list, sep = "\n")

files_list = get_files_list(year,version,type)
samples_type = [os.path.basename(f).replace('_rootFileList.txt','') for f in files_list]



# os.chdir(work_dir)

# options for the macro
#   -w [ --workflow ] arg (=0)   Workflow number
#   --workflow_title arg         Workflow title
#   --mc                         Run on Monte Carlo
#   --data                       Run on data

macro_opts  = f'--workflow={workflow} '
macro_opts += f'--workflow_title=\"{workflow_title}\" '
macro_opts += f'--{type}'
macro_opts  = f'\'{macro_opts}\''


# options for naf_submit
#   --exe EXE, -e EXE              Executable (REQUIRED)
#   --config CONFIG, -c CONFIG     Configuration file (REQUIRED)
#   --ntuples NTUPLES, -n NTUPLES  List of ntuples file
#   --output OUTPUT, -o OUTPUT     Name of output root file
#   --nfiles NFILES, -x NFILES     Number of ntuple files per job
#   --json JSON, -j JSON           JSON file with certified data
#   --label LABEL, -l LABEL        user label for the submission
#   --opts OPTS                    additional options for the macro execution (job.sh)
#   --events EVENTS_MAX            override eventsMax in the config file (default = -1)

# QCD samples for each trigger
qcd_sample = {}
qcd_sample['Jet40'] = ['QCD_HT100to200','QCD_HT200to300','QCD_HT300to500']
qcd_sample['Jet100'] = ['QCD_HT200to300','QCD_HT300to500','QCD_HT500to700','QCD_HT700to1000']
qcd_sample['Jet200'] = ['QCD_HT500to700','QCD_HT700to1000','QCD_HT1000to1500','QCD_HT1500to2000','QCD_HT2000toInf']
qcd_sample['Jet350'] = ['QCD_HT1500to2000','QCD_HT2000toInf']

naf_opts = ''
for alias,cfg in config_files.items():
    for i,file_list in enumerate(files_list):
        # if i != 3: continue
        # if type == 'mc' and samples_type[i] not in qcd_sample[alias]: continue
        output_root = f'histograms_{year}_{version}_{type}_{samples_type[i]}_{alias}.root'
        naf_opts  = f' --exe=TagAndProbeOnlineBtag'
        naf_opts += f' --config={cfg}'
        naf_opts += f' --ntuples={file_list}'
        naf_opts += f' --nfiles=5'
        naf_opts += f' --label=wf{workflow}_{samples_type[i]}'
        naf_opts += f' --events=-1'
        naf_opts += f' --output={output_root}'
        # naf_opts += f' --json={json}'  # TODO: investigate the proble with json; hard-code in the template
        naf_opts += f' --opts={macro_opts}'

        naf_command = f'cd {work_dir} ; naf_submit.py {naf_opts}'
        
        print(naf_command)
        os.system(naf_command)
        print()




