#!/usr/bin/env python3
import os
import shutil
import argparse
import ROOT

# This script is used to add histograms from DQM .root files.
# It will search for all .root files in the source directory, extract the histograms for a given trigger,
# and combine them into a single .root file.
# The histograms are expected to be in the format:
# DQMData/Run <run_number>/HLT/Run summary/SUS/MssmHbb/control/btag/<trigger_name>/<histogram_name>
# The script will create a new .root file with the combined histograms.
# The script can be run from the command line with the following arguments:
# - source_dir: The directory containing the .root files to be processed
# - --target_dir: Optional argument to specify the target directory for the output file
# - --trigger_name: The name of the trigger to be processed, which will be used to locate the histograms in the .root files.    
# example usage: ./add_histograms.py /eos/home-r/rwalsh/ww3/battle/4/trigger_studies_2025/efficiencies_dqm/data --target_dir $PWD/2025C --trigger_name HLT_DoublePFJets40_PNetBTag_0p11


def get_histograms(input_file, run_number, trigger_name, histogram_list):
    histograms = {}
    
    for h_name in histogram_list:
        histograms[h_name] = input_file.Get(f"DQMData/Run {run_number}/HLT/Run summary/SUS/MssmHbb/control/btag/{trigger_name}/{h_name}")

    return histograms


def main():
    
    histograms_list = []
    histograms_list.append("pt_probe")
    histograms_list.append("pt_probe_match")
    histograms_list.append("discr_offline_btag_jet1")
    histograms_list.append("discr_offline_btag_jet2")
    
    parser = argparse.ArgumentParser(description="Add histograms from DQM .root files.")
    parser.add_argument("source_dir", help="Source directory containing .root files")
    parser.add_argument("--target_dir", help="Optional target directory for copied files")
    parser.add_argument("--trigger_name", help="Root directory containing the histograms.")

    args = parser.parse_args()
    
    combined = {}
    for hist_name in histograms_list:
        combined[hist_name] = None
    
    # combined_pt_probe = None
    # combined_pt_probe_match = None
    # combined_offline_btag_jet1 = None
    # combined_offline_btag_jet2 = None
    
    count_files = 0
    
    for root_path, dirs, files in os.walk(args.source_dir):
        for filename in files:
            if filename.endswith(".root") and "PromptReco" in filename:
                count_files += 1
                full_filename = os.path.join(root_path, filename)
                run_number = filename.split("_")[2].replace("R000","")
                input_file = ROOT.TFile.Open(full_filename, "READ")
                # pt_probe, pt_probe_match, discr_offline_btag_jet1, discr_offline_btag_jet2 = get_histograms(input_file, run_number, args.trigger_name)
                histograms = get_histograms(input_file, run_number, args.trigger_name, histograms_list)
                
                print(f"{count_files:03d},{run_number},{input_file}")
                
                for hist_name, histogram in histograms.items():
                    if combined[hist_name] is None:
                        combined[hist_name] = histogram.Clone()
                        combined[hist_name].SetDirectory(0)
                    else:
                        combined[hist_name].Add(histogram)
                
                # if combined_pt_probe is None:
                #     combined_pt_probe = pt_probe.Clone()
                #     combined_pt_probe.SetDirectory(0)
                #     combined_pt_probe_match = pt_probe_match.Clone()
                #     combined_pt_probe_match.SetDirectory(0)
                #     combined_offline_btag_jet1 = discr_offline_btag_jet1.Clone()
                #     combined_offline_btag_jet1.SetDirectory(0)
                #     combined_offline_btag_jet2 = discr_offline_btag_jet2.Clone()
                #     combined_offline_btag_jet2.SetDirectory(0)
                # else:
                #     combined_pt_probe.Add(pt_probe)
                #     combined_pt_probe_match.Add(pt_probe_match)
                #     combined_offline_btag_jet1.Add(discr_offline_btag_jet1)
                #     combined_offline_btag_jet2.Add(discr_offline_btag_jet2)
                    

    if args.target_dir:
        os.makedirs(args.target_dir, exist_ok=True)
        output_file_path = os.path.join(args.target_dir, f"btag_probe_{args.trigger_name}.root")
    else:
        output_file_path = f"btag_probe_{args.trigger_name}.root"
    
    output_file = ROOT.TFile(output_file_path, "recreate")
    for hist_name in histograms_list:
        combined[hist_name].Write()
    
    # combined_pt_probe.Write()
    # combined_pt_probe_match.Write()
    # combined_offline_btag_jet1.Write()
    # combined_offline_btag_jet2.Write()
    output_file.Close()
    
if __name__ == "__main__":
    main()
