#python prepare_for_multibeam_folds.py --cands /b/PROCESSING/12_CANDYJAR/3HM_ACCELSEARCH_FULL/TOP_300/classifications/T1_CANDS.csv --out_prefix 3HM_T1 --beam_root /b/PROCESSING/02_FILTERBANKS/J0514-4002A/2024-05-19-15:50:23/ --overlap_csv /b/PROCESSING/01_BEAMFORMED/J0514-4002A/2024-05-19-15:50:23/swdelays_J0514-4002A_1716133768_to_1716140969_f63918_overlapping_beams.csv  --output_root /bscratch/CAND_RECHECKS_6
import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import text
import uuid
import xml.etree.ElementTree as ET
import os
import argparse
import configparser
import ast
import glob
from itertools import chain
# Replace with your actual connection string
speed_of_light = 299792458.0

database_ini = os.environ.get('HOME') + '/.compactdb_secrets.ini'

config = configparser.ConfigParser()
config.read(database_ini)

DB_USER = config['compactdb_ro']['DB_USER']
DB_PASSWORD = config['compactdb_ro']['DB_PASSWORD']
DB_HOST = config['compactdb_ro']['DB_HOST']
DB_PORT = config['compactdb_ro']['DB_PORT']
DB_NAME = config['compactdb_ro']['DB_NAME']

DATABASE_URI = (f"mysql+pymysql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}")
engine = sqlalchemy.create_engine(DATABASE_URI)



class XMLValues(object):
    def __init__(self, fil_file, tsamp, fft_size, coherent_dm, tstart, segment_pepoch):
        self.fil_file = fil_file
        self.tsamp = tsamp
        self.fft_size = fft_size
        self.coherent_dm = coherent_dm
        self.tstart = tstart
        self.segment_pepoch = segment_pepoch

def period_correction_for_prepfold(p0,pdot,tsamp,fft_size):
    return p0 - pdot*fft_size*tsamp/2.

def f0_correction_for_prepfold(f0,fdot,tsamp,fft_size):
    return f0 - fdot*float(fft_size)*tsamp/2.


def calculate_spin(f=None, fdot=None, p=None, pdot=None):
        # calculate p and pdot from f and fdot
        if f is not None and fdot is not None:
            p = 1 / f
            pdot = -fdot / (f**2)
        # calculate f and fdot from p and pdot
        elif p is not None and pdot is not None:
            f = 1 / p
            fdot = -pdot / (p**2)
        else:
            raise ValueError("Either (f, fdot) or (p, pdot) must be provided")
        
        a = pdot / p * speed_of_light
            
        return f, fdot, p, pdot, a

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare for multibeam folds - generate prepfold commands to fold same beams, psrfold commands to fold all overlapping and neighbouring beams")
    parser.add_argument(
        "--cands",
        type=str,
        required=True,
        help="Path to cands file"
    )
    parser.add_argument(
        "--out_prefix",
        type=str,
        required=True,
        help="Prefix for output files"
    )

    parser.add_argument(
        "--beam_root",
        type=str,
        required=True,
        help="Root directory for beams",
    )

    parser.add_argument(
        "--output_root",
        type=str,
        required=True,
        help="Root directory for folds",

    )

    parser.add_argument(
        "--overlap_csv",
        type=str,
        required=True,
        help="Path to overlap csv file"
    )

    parser.add_argument(
        "--use_search",
        action='store_true',
        required=False,
        default=False,
        help="Use search candidate details instead of optimised fold values"
    )

    return parser.parse_args()



if __name__ == "__main__":
    args = parse_args()
    os.makedirs(args.output_root, exist_ok=True)
    t1_df = pd.read_csv(args.cands)
    fold_ids = t1_df['fold_candidate_id_bin']

    print(f"Successfully read {len(t1_df.index)} T1 candidates")

    fold_ids_hex = [f"\"{uuid.UUID(x).hex}\"" for x in fold_ids]
    fold_id_list_str = ", ".join(fold_ids_hex)
    fold_candidate_df = pd.read_sql_query(f"""
        SELECT 
        HEX(f.id) AS foldcand_hex,
        f.spin_period, f.pdot, f.dm, f.fold_snr, f.filepath, f.filepath, HEX(f.dp_id) AS dp_hex, HEX(f.search_candidate_id) AS sc_hex
                                        
        FROM fold_candidate f
        WHERE HEX(f.id) IN ({fold_id_list_str})
    """, engine)

    print("Obtained relevant fold candidate information from database")



    foldcand_xml_map_df = pd.read_sql_query(text(f"""
        SELECT 
        HEX(f.id)           AS foldcand_hex,
        HEX(dp_child.id)    AS child_dp_hex,
        HEX(pdi.dp_id)      AS xml_dp_hex
        FROM fold_candidate f
        JOIN data_product dp_child 
            ON dp_child.id = f.dp_id
        JOIN processing_dp_inputs pdi
            ON pdi.task_id  = dp_child.created_by_task_id
            AND pdi.run_name= dp_child.created_by_run_name
            AND pdi.session_id= dp_child.created_by_session_id
        WHERE HEX(f.id) IN ({fold_id_list_str})
    """), engine)
    

    unique_xml_dp_hex = foldcand_xml_map_df['xml_dp_hex'].unique()
    if len(unique_xml_dp_hex) == 0:
        # Edge case: no parent XML. Possibly raise an exception or just proceed.
        # We'll proceed with an empty result if that happens.
        xml_dps_df = pd.DataFrame()
    else:
        xml_dp_list_str = ", ".join(f"\"{x}\"" for x in unique_xml_dp_hex)
        xml_dps_df = pd.read_sql_query(f"""
            SELECT
            HEX(dp.id) AS xml_dp_hex,
            dp.filepath AS xml_filepath,
            dp.filename AS xml_filename
            FROM data_product dp
            WHERE HEX(dp.id) IN ({xml_dp_list_str})
        """, engine)


    foldcand_xml_merged = pd.merge(
        foldcand_xml_map_df,
        xml_dps_df,
        on='xml_dp_hex',
        how='left'
    )

    print("Obtained corresponding XML names from database")

    foldcand_to_fil = {}

    for row in foldcand_xml_merged.itertuples(index=False):
        foldcand_hex = row.foldcand_hex
        xml_fullpath = os.path.join(row.xml_filepath, row.xml_filename)
        
        if not os.path.isfile(xml_fullpath):
            # Handle the possibility that file is missing
            foldcand_to_fil[foldcand_hex] = None
            continue
        
        tree = ET.parse(xml_fullpath)
        root = tree.getroot()

        fil_file = root.find('search_parameters/infilename')
        tsamp = root.find('header_parameters/tsamp')
        fft_size = root.find('search_parameters/size')
        coherent_dm = root.find('search_parameters/cdm')
        tstart = root.find('header_parameters/tstart')
        segment_pepoch = root.find('segment_parameters/segment_pepoch')

        if tsamp is not None and fft_size is not None and tstart is not None and segment_pepoch is not None and fil_file is not None and coherent_dm is not None:
            foldcand_to_fil[foldcand_hex] = XMLValues(fil_file.text, float(tsamp.text), int(fft_size.text), float(coherent_dm.text), float(tstart.text), float(segment_pepoch.text))
        else:
            print(f"Missing required field in XML: {xml_fullpath}, values are tsamp: {tsamp}, fft_size: {fft_size}, tstart: {tstart}, segment_pepoch: {segment_pepoch}, fil_file: {fil_file}, coherent_dm: {coherent_dm}")
            continue



    fold_candidate_df['fil_file'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].fil_file if x in foldcand_to_fil else None)
    fold_candidate_df['tsamp'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].tsamp if x in foldcand_to_fil else None)
    fold_candidate_df['fft_size'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].fft_size if x in foldcand_to_fil else None)
    fold_candidate_df['coherent_dm'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].coherent_dm if x in foldcand_to_fil else None)
    fold_candidate_df['tstart'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].tstart if x in foldcand_to_fil else None)
    fold_candidate_df['segment_pepoch'] = fold_candidate_df['foldcand_hex'].apply(lambda x: foldcand_to_fil[x].segment_pepoch if x in foldcand_to_fil else None)

    t1_df['fold_ids_hex'] = t1_df['fold_candidate_id_bin'].apply(lambda x: uuid.UUID(x).hex.upper())
    fold_candidate_df = fold_candidate_df.merge(t1_df[['fold_ids_hex', 'f0_usr', 'f1_user', 'dm_user', 'beam_name']], left_on='foldcand_hex', right_on='fold_ids_hex', how='inner')

    #rename user values to search
    fold_candidate_df.rename(columns={'f0_usr': 'f0_search', 'f1_user': 'f1_search', 'dm_user': 'dm_search'}, inplace=True)



    fold_candidate_df["corrected_spin_period"]  = period_correction_for_prepfold(fold_candidate_df['spin_period'], fold_candidate_df['pdot'], fold_candidate_df['tsamp'], fold_candidate_df['fft_size'])


    fold_candidate_df.rename(columns={'beam_name': 'beam'}, inplace=True)

    print("Successfully merged fold candidate information with XML data")

    if args.use_search:

        fold_candidate_df["corrected_f0_search"]  = fold_candidate_df.apply(
            lambda row: f0_correction_for_prepfold(row['f0_search'], row['f1_search'], row['tsamp'], row['fft_size']),
            axis=1
        )


        fold_candidate_df["prepfold_cmd"] = fold_candidate_df.apply(
            lambda row: (
                f"prepfold -topo -ncpus 2 -fixchi -noxwin"
                f" -dm {row['dm_search']} "
                f" -n 64 -nsub 64 "
                f" -f {row['corrected_f0_search']}"
                f" -fd {row['f1_search']} "
            ),
            axis=1
        )
    else:
        fold_candidate_df["prepfold_cmd"] = fold_candidate_df.apply(
            lambda row: (
                f"prepfold -topo -ncpus 2 -fixchi -noxwin"
                f" -dm {row['dm']} "
                f" -n 64 -nsub 64 "
                f" -p {row['corrected_spin_period']}"
                f" -pd {row['pdot']} "
            ),
            axis=1
        )


    value_suffix = 'search-vals' if args.use_search else 'fold-vals'

    csv_df = fold_candidate_df 
    prepfold_csv_name = args.output_root + "/" + args.out_prefix +'_prepfold_cmds_'  + value_suffix  + '.csv'
    csv_df.to_csv(prepfold_csv_name, index=False)

    print(f"Successfully wrote prepfold commands to {prepfold_csv_name}")


    overlap_df = pd.read_csv(args.overlap_csv)
    overlap_df = overlap_df.rename( columns={'name': 'beam'})

    print(f"Successfully read overlap csv file with {len(overlap_df.index)} rows")

    overlap_df['all_related_beams'] = overlap_df.apply(
        lambda row: ([row['beam']]) + ast.literal_eval(row['overlapping_beams'] or []) + ast.literal_eval(row['neighbouring_beams'] or []) ,
        axis=1
    )

    overlap_df.drop(columns=['overlapping_beams','neighbouring_beams' ])

    total_beams = np.unique(np.array(list(chain.from_iterable(overlap_df['all_related_beams']))))
    merged_df = fold_candidate_df.merge(
        overlap_df[['beam', 'all_related_beams']],
        on='beam',
        how='left'
    )

    total_beams = np.unique(np.array(list(chain.from_iterable(overlap_df['all_related_beams']))))
    total_dms = np.unique(merged_df['coherent_dm'])
    folds_per_beam_per_dm = {}
    for beam in total_beams:
        for dm in total_dms:
            #shortlist merged_df if coherent_dm ==dm and beam in all_related_beams list
            shortlist = merged_df[(merged_df['coherent_dm']==dm) & (merged_df['all_related_beams'].apply(lambda x: beam in x))]
            #print(f"beam {beam} and dm {dm} has {len(shortlist.index)} related candidates")
            if len(shortlist.index) > 0:
                folds_per_beam_per_dm[(beam,dm)] = shortlist

    print(f"Successfully grouped fold candidates by beam and DM")

    #make psrfold commands
    beam_root=args.beam_root
    fold_root = args.output_root
    os.makedirs(fold_root, exist_ok=True)
    psrfold_df = pd.DataFrame(columns=['beam', 'dm', 'psrfold_cmd'])

    idx=0
    for key, value in folds_per_beam_per_dm.items():
        beam, dm = key
        candfile = os.path.join(fold_root, f"{beam}_{dm:.02f}_{value_suffix}.candfile")
        with open(candfile, 'w') as cand_file_writer:
            cand_file_writer.write("#id DM accel F0 F1 F2 S/N\n")
            if args.use_search:
                for i, row in value.iterrows():
                    f, fdot, p, pdot, a = calculate_spin(f=row['f0_search'], fdot=row['f1_search'])
                    cand_file_writer.write(f"{i} {row['dm_search']} {a} {f} {fdot} 0 {row['fold_snr']}\n")
            else:
                for i, row in value.iterrows():
                    f, fdot, p, pdot, a = calculate_spin(p=row['corrected_spin_period'], pdot=row['pdot'])
                    cand_file_writer.write(f"{i} {row['dm']} {a} {f} {fdot} 0 {row['fold_snr']}\n")

        fil_files = glob.glob(f"{beam_root}/*/{beam}/*idm_{dm:09.3f}_{beam}*.fil")
        psrfold_cmd = f"psrfold_fil --plotx -v --render --candfile {candfile} \
        -z zap 552.0 553.0 -z zap 926.0 927.0 -z zap 934.0 952.0 -z zap 1036.0 1037.0 -z zap 1062.0 1063.0 -z zap 1079.0 1080.0    \
        --template /home/pulsarx/software/PulsarX/include/template/meerkat_fold_UHF.template --clfd 2.0 --rfi zdot --pepoch  {row['segment_pepoch']} -n 64 -b 128 -f {' '.join(fil_files)}"
        psrfold_df.loc[idx] = [beam, dm, psrfold_cmd]
        idx+=1

    #sort by beam
    psrfold_df = psrfold_df.sort_values(by=['beam'])
    psrfold_csv_name = args.output_root + "/" + args.out_prefix +'_psrfold_cmds.csv'
    psrfold_df.to_csv(psrfold_csv_name, index=False)

    print(f"Successfully wrote psrfold commands to {psrfold_csv_name}")

    print("All done!")

    


        
                



