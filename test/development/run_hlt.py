#!/bin/env python3

import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description='Read a text file containing a list of data files.')
    parser.add_argument('--campaign', type=str, help='the campaign parameter')
    parser.add_argument('--mass', type=str, help='the mass parameter')
    parser.add_argument('--num_files', type=int, help='the number of data files to read')
    parser.add_argument('--menu_path', type=str, default='/users/rwalsh/run3/mssmhbb/V135', help='menu path')    
    args = parser.parse_args()
    
    filename = f'samples/{args.campaign}_M-{args.mass}.dat'
    
    if not os.path.exists('nfs'):
        print("Please create an 'nfs' directory and link it here.")
        sys.exit()

    prefix = "file:/pnfs/desy.de/cms/tier2"
    with open(filename, 'r') as f:
        data_files = [f'{prefix}{line.strip()}' for line in f.readlines()[:args.num_files]]
    files_list = ','.join(data_files)
    
    menu_version = args.menu_path.split('/')[-1].lower()
    
    run_dir = f'nfs/{args.campaign}/{menu_version}'
    os.makedirs(run_dir,exist_ok=True)
    os.chdir(run_dir)    
    
    hlt_get_config = f'''hltGetConfiguration {args.menu_path}   \\
--setup /dev/CMSSW_13_3_0/GRun/V14  \\
--globaltag auto:phase1_2023_realistic    \\
--mc    \\
--unprescale    \\
--output minimal    \\
--max-events -1    \\
--input {files_list}   \\
--eras Run3 \\
--l1-emulator FullMC \\
--l1 L1Menu_Collisions2023_v1_3_0_xml  \\
--dbtunnelport 10121   \\
--dbtunnel > hltMC_mssmhbb_run3_m{args.mass}.py'''

    nice_run = f'( nice cmsRun hltMC_mssmhbb_run3_m{args.mass}.py &> hltMC_mssmhbb_run3_m{args.mass}.log & )'
    
    os.system(f'{hlt_get_config}; {nice_run}')
    # subprocess.run([hlt_get_config + ' && ' + nice_run], shell=True)

if __name__ == '__main__':
    main()
