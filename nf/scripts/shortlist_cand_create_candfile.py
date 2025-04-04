import xml.etree.ElementTree as ET
import sys, os, subprocess
import argparse
import pandas as pd
import numpy as np
import logging
import time
import shlex
import threading
from multiprocessing import Pool, cpu_count
import re
import json

def generate_pulsarX_cand_file(out_dir, global_filtered_df):
    cand_dms = global_filtered_df['dm'].values
    cand_accs = global_filtered_df['acc'].values
    cand_period = global_filtered_df['period'].values
    cand_mod_frequencies = 1 / cand_period
    cand_snrs = global_filtered_df['snr'].values
    cand_file_path = os.path.join(out_dir, 'output.candfile')
    with open(cand_file_path, 'w') as f:
        f.write("#id DM accel F0 F1 F2 S/N\n")
        for i in range(len(cand_mod_frequencies)):
            f.write("%d %f %f %f 0 0 %f\n" % (i, cand_dms[i], cand_accs[i], cand_mod_frequencies[i], cand_snrs[i]))

    return cand_file_path

def get_shortlisted_candidate_df(input_cand_df, candidates, ptol, dm_tol, ncands, birdies):
    ignored_entries = [
        'candidate', 'opt_period', 'folded_snr', 'byte_offset', 'is_adjacent',
        'is_physical', 'ddm_count_ratio', 'ddm_snr_ratio'
    ]
    rows = []
    for candidate in candidates:
        cand_dict = {}
        for cand_entry in candidate.iter():
            if cand_entry.tag not in ignored_entries:
                cand_dict[cand_entry.tag] = cand_entry.text
        cand_dict['cand_id_in_file'] = candidate.attrib.get("id")
        rows.append(cand_dict)

    xml_df = pd.DataFrame(rows)
    xml_df = xml_df.astype({"snr": float, "dm": float, "period": float, "nh": int, "acc": float, "nassoc": int, "cand_id_in_file": int})

    if birdies:
        tol=1e-6
        birdies = [float(birdie) for birdie in birdies.split(',')]
        # Filter out candidates with frequencies close to the birdies
        for birdie in birdies:
            xml_df = xml_df[~xml_df['period'].astype(float).between(birdie-tol, birdie+tol)]
            #xml_df = xml_df[~np.isclose(xml_df['period'].astype(float), birdie, atol=tol)]

    global_filtered_df = pd.DataFrame(columns=xml_df.columns)

    for _,row in input_cand_df.iterrows():
        F0=float(row['F0'])
        period = 1/F0
        DM=float(row['DM'])
        min_period = period - (period * ptol)
        max_period = period + (period * ptol)
    
        min_dm = DM - dm_tol
        max_dm = DM + dm_tol

        # Filter the xml_df based on the period and DM
        filtered_df = xml_df[
            (xml_df['period'].astype(float) >= min_period) &
            (xml_df['period'].astype(float) <= max_period) &
            (xml_df['dm'].astype(float) >= min_dm) &
            (xml_df['dm'].astype(float) <= max_dm)
        ]
        # Add the filtered dataframe to the global filtered dataframe
        global_filtered_df = pd.concat([global_filtered_df, filtered_df], ignore_index=True)
        
    global_filtered_df = global_filtered_df.drop_duplicates(subset=['cand_id_in_file'])

    # get remaining candidates from xml_df sorted by snr to reach max number of candidates to fold
    remaining_candidates = xml_df[
        ~xml_df['cand_id_in_file'].isin(global_filtered_df['cand_id_in_file'])
    ].sort_values(by='snr', ascending=False)

    print(f"{len(global_filtered_df)} candidates shortlisted from xml file")
    print(f"{len(remaining_candidates)} candidates remaining in xml file")


    # get the number of candidates to add from remaining candidates
    num_to_add = ncands - len(global_filtered_df)
    if num_to_add > 0:
        remaining_candidates = remaining_candidates.head(num_to_add)
        global_filtered_df = pd.concat([global_filtered_df, remaining_candidates], ignore_index=True)
    # sort the final dataframe by snr
    global_filtered_df = global_filtered_df.sort_values(by='snr', ascending=False)
    return global_filtered_df

def write_pulsarx_commands(
    out_dir, pepoch, input_filenames, source_name_prefix,
    nbins_high, nbins_low, subint_length, nsubband,
    beam_name, pulsarx_threads, template, clfd_q_value,
    rfi_filter, start_fraction=None, end_fraction=None,
     coherent_dm=0.0):
    beam_tag = "-i {}".format(int(beam_name.replace("cfbf", "").strip()))
    if rfi_filter:
        additional_flags = f"--rfi {rfi_filter}"
    else:
        additional_flags = ""
    nbins_string = "-b {} --nbinplan 0.01 {}".format(nbins_low, nbins_high)
    commands= []
    for pulsarx_predictor, output_dir in zip(['input.candfile', 'output.candfile'], ['direct_fold', 'search_fold']):

        output_rootname = output_dir + "/" + beam_name

        script = (
            "psrfold_fil2 -v --render --output_width --cdm {} -t {} --candfile {} -n {} {} {}  --template {} "
            "--clfd {} -L {} -f {} -o {} --srcname {} --pepoch {} --frac {} {} {}"
        ).format(
            coherent_dm,
            pulsarx_threads,
            pulsarx_predictor,
            nsubband,
            nbins_string,
            beam_tag,
            template,
            clfd_q_value,
            subint_length,
            input_filenames,
            output_rootname,
            source_name_prefix,
            pepoch,
            start_fraction,
            end_fraction,
            additional_flags
        )
        commands.append(script)

    
    with open(os.path.join(out_dir, "pulsarx_cmds.txt"), "w") as f:
        for cmd in commands:
            f.write(cmd + "\n")


def main():
    parser = argparse.ArgumentParser(description="Shortlist and Fold")
    parser.add_argument(
        "--xml", type=str, required=True, help="Input XML file"
    )
    parser.add_argument(
        "--input_candfile",
        type=str,
        required=True,
        help="Input candidate file")
    parser.add_argument(
        "--ptol",
        type=float,
        default=0.01,
        help="Fold tolerance fraction of spin period (default: 0.01)")
    parser.add_argument(
        "--dm_tol",
        type=float,
        default=10,
        help="DM tolerance in absolute units (default: 10)")
    
    parser.add_argument(
        "--ncands",
        type=int,
        default=1000,
        help="Number of candidates to shortlist (default: 1000)")
    parser.add_argument(
        "--out_dir",
        type=str,
        default=".",
        help="Output directory (default: current directory)")
    
    parser.add_argument(
        "--birdies",
        type=str, 
        default="0.004615,0.179999,0.18000,0.004614,0.004616,0.001513",
        help="Comma-separated list of birdies to ignore (default: empty string)"

    )
    parser.add_argument(
        "--template_dir",
        type=str,
        required=True,
        help="Directory containing the template files"
    )
    parser.add_argument('-u', '--nbins_high', help='Upper profile bin limit for slow-spinning pulsars',
                        type=int, default=64)
    parser.add_argument('-l', '--nbins_low', help='Lower profile bin limit for fast-spinning pulsars',
                        type=int, default=32)
    parser.add_argument('-sub', '--subint_length', help='Subint length (s). Default is tobs/64',
                        type=int, default=None)
    parser.add_argument('-nsub', '--nsubband', help='Number of subbands',
                        type=int, default=64)
    parser.add_argument('-clfd', '--clfd_q_value', help='CLFD Q value',
                        type=float, default=2.0)
    parser.add_argument('-rfi', '--rfi_filter', help='RFI filter value',
                        type=str, default=None)
    parser.add_argument('-b', '--beam_name', help='Beam name string',
                        type=str, default='cfbf00000')
    parser.add_argument('-threads', '--pulsarx_threads', help='Number of threads to be used for pulsarx',
                        type=int, default=1)
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose (DEBUG) logging')
    parser.add_argument('--extra_args', type=str, default=None,
                        help='Extra arguments to pass to psrfold_fil.')
    parser.add_argument('--cdm', type=float, default=0.0,
                        help='Coherent DM to use for folding. Default is 0.0.')
    args = parser.parse_args()

    input_candfile = args.input_candfile

    input_candfile = os.path.realpath(os.path.abspath(input_candfile))
    #copy to output dir
    os.system(f"cp {input_candfile} {args.out_dir}/input.candfile")
    
    #candfile is of the space separated format, get dataframe from it
    input_cand_df = pd.read_csv(input_candfile, delim_whitespace=True)

    print(input_cand_df.columns)


    xml_file = args.xml
    tree = ET.parse(xml_file)
    root = tree.getroot()
    header_params = root[1]    
    search_params = root[2]
    segment_params = root[3]
    candidates = root[7]

    filterbank_file = str(search_params.find("infilename").text)
    filterbank_file = os.path.realpath(os.path.abspath(filterbank_file))
    if not os.path.exists(filterbank_file):
        print("Filterbank file not found: ", filterbank_file)
        sys.exit(1)
    os.makedirs(args.out_dir, exist_ok=True)

    segment_start_sample = int(segment_params.find('segment_start_sample').text)
    segment_nsamples = int(segment_params.find('segment_nsamples').text)
    segment_pepoch = float(segment_params.find('segment_pepoch').text)

    tsamp = float(header_params.find("tsamp").text)
    fch1 = float(header_params.find("fch1").text)
    source_name_prefix = str(header_params.find("source_name").text).strip()

    template = None
    if fch1 < 1100:
        template = os.path.join(args.template_dir,  "meerkat_fold_UHF.template")    
    elif fch1 >= 1700 and fch1 <= 1712:
        template = os.path.join(args.template_dir, "meerkat_fold_L.template")
    elif fch1 > 1712:
        template = os.path.join(args.template_dir, "meerkat_fold_S.template")
    else:
        print("Frequency not in range")
        sys.exit(1)
    total_nsamples = int(root.find("header_parameters/nsamples").text)
    user_start_fraction = round(segment_start_sample / total_nsamples, 3)
    user_end_fraction = round((segment_start_sample + segment_nsamples) / total_nsamples, 3)
    
    effective_tobs = tsamp * segment_nsamples

    print("Effective tobs: ", effective_tobs)

    if args.subint_length is None:
        subint_length = int(effective_tobs / 64)
    else:
        subint_length = args.subint_length

    global_filtered_df = get_shortlisted_candidate_df(input_cand_df, candidates, args.ptol, args.dm_tol, args.ncands, args.birdies)
    generate_pulsarX_cand_file(args.out_dir, global_filtered_df)
    write_pulsarx_commands(
        args.out_dir,
        segment_pepoch,
        filterbank_file,
        source_name_prefix,
        args.nbins_high,
        args.nbins_low,
        subint_length,
        args.nsubband,
        args.beam_name,
        args.pulsarx_threads,
        template,
        args.clfd_q_value,
        args.rfi_filter,
        user_start_fraction,
        user_end_fraction,
        coherent_dm=args.cdm
    )

    print("Shortlisted candidates:")
    print(global_filtered_df[['cand_id_in_file', 'snr', 'dm', 'period', 'acc']])
    print("Number of shortlisted candidates: ", len(global_filtered_df))
    print("All done")



if __name__ == "__main__":
    main()