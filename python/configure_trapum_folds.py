import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import coordinates as coord
from astropy import wcs
from matplotlib.patches import Ellipse
import math
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import json
import argparse
from matplotlib import rcParams
params = {"figure.figsize": (12, 9),
          "font.size": 12,
          "font.family": 'serif',
          "font.weight": "normal",
          "xtick.major.size": 9,
          "xtick.minor.size": 4,
          "ytick.major.size": 9,
          "ytick.minor.size": 4,
          "xtick.major.width": 1,
          "xtick.minor.width": 1,
          "ytick.major.width": 1,
          "ytick.minor.width": 1,
          "xtick.major.pad": 8,
          "xtick.minor.pad": 8,
          "ytick.major.pad": 8,
          "ytick.minor.pad": 8,
          "lines.linewidth": 1.2,
          "lines.markersize": 10,
          "axes.linewidth": 2.5,
          "legend.loc": "upper right",
          "text.usetex": False,
          "xtick.labelsize": 12,
          "ytick.labelsize": 12,
          "xtick.direction": "in",
          "ytick.direction": "in",
          "xtick.bottom": True,
          "xtick.top": True,
          "ytick.left": True,
          "ytick.right": True,
          }

rcParams.update(params)


def convert_equatorial_coordinate_to_pixel(equatorial_coordinates, bore_sight, utc_time):

    """
    https://docs.astropy.org/en/stable/wcs/index.html#using-astropy-wcs
    """
    step = 1/10000000000.

    wcs_properties = wcs.WCS(naxis=2)
    wcs_properties.wcs.crpix = [0, 0]
    wcs_properties.wcs.cdelt = [-step, step]
    wcs_properties.wcs.crval = [bore_sight.ra.deg,bore_sight.dec.deg]
    wcs_properties.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    scaled_pixel_coordinats = wcs_properties.wcs_world2pix([[equatorial_coordinates.ra.deg,equatorial_coordinates.dec.deg]], 0)
    pixel_coordinates = scaled_pixel_coordinats * step

    return pixel_coordinates

def get_pixel_coherent_beam_coordinates(beam_coords, boresight_coords, utc_time):
    """
    Convert coherent beam equatorial coordinates to pixel coordinates
    """
    pixel_beam_ras=[]
    pixel_beam_decs=[]
    # Convert equatorial beam coordinates to pixel coordinates
    for beam_coord in beam_coords:
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(beam_coord, boresight_coords, utc_time)
        pixel_beam_ras.append(boresight_coords.ra.deg + pixel_coordinates[0][0])
        pixel_beam_decs.append(boresight_coords.dec.deg + pixel_coordinates[0][1])

    return pixel_beam_ras, pixel_beam_decs

def ellipse_parametric(t, a, b, x0, y0, theta):
        cos_t = np.cos(t)
        sin_t = np.sin(t)
        x = x0 + a * cos_t * np.cos(theta) - b * sin_t * np.sin(theta)
        y = y0 + a * cos_t * np.sin(theta) + b * sin_t * np.cos(theta)
        return x, y

def position_within_beam(beam, pixel_ra, pixel_dec):
    x = pixel_ra
    y = pixel_dec
    x0, y0, a, b, theta = beam.pixel_ra, beam.pixel_dec, beam.x, beam.y, beam.angle
    cos_theta = np.cos(-theta)
    sin_theta = np.sin(-theta)
    xr = cos_theta * (x - x0) - sin_theta * (y - y0)
    yr = sin_theta * (x - x0) + cos_theta * (y - y0)
    return (xr**2 / a**2) + (yr**2 / b**2) <= 1

def check_containment( beam1, beam2):
    # Check if the center of ellipse1 is inside ellipse2
    return position_within_beam(beam2, beam1.pixel_ra, beam1.pixel_dec) or  position_within_beam(beam1, beam2.pixel_ra, beam2.pixel_dec)

def discrete_overlap( beam1, beam2, num_points=100):
    if check_containment(beam1, beam2):
        return True
    # Check points on the perimeters of both ellipses
    t_values = np.linspace(0, 2 * np.pi, num_points)
    #psf_x, psf_y -> semi-major axis (a), semi-minor axis (b) of the ellipse
    x1, y1 = ellipse_parametric(t_values, beam1.x, beam1.y,beam1.pixel_ra, beam1.pixel_dec, beam1.angle)
    
    # Check if any points on the perimeter of ellipse1 lie inside ellipse2
    for x, y in zip(x1, y1):
        if position_within_beam(beam2, x, y):
            return True    
    return False



def parse_meta_file(meta_file):
    with open(meta_file) as f:
        meta_file_json = json.load(f)
        return meta_file_json
    
class Beam(object):
    def __init__(self, name, utc, pid, freq, target, ra, dec,pixel_ra, pixel_dec, x, y, angle):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.pixel_ra = pixel_ra
        self.pixel_dec = pixel_dec
        self.x = x
        self.y = y
        self.angle = angle
        self.utc = utc
        self.overlapping_beams = []
        self.filterbank_list = []
        self.pid = pid
        self.freq = freq
        self.target = target
    
    def add_overlapping_beams(self, beams):
        self.overlapping_beams.extend(beams)
    
    def add_filterbanks(self, filterbanks):
        self.filterbank_list.extend(filterbanks)
    
    def __str__(self):
        return f"Beam {self.name} "
    
    def __repr__(self):
        return f"Beam {self.name} "
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)
        

#python configure_trapum_folds.py --meta1 "/b/PROCESSING/01_BEAMFORMED/J0514-4002A/2024-05-19-15:50:23/3HM.meta" --meta2 "/bscratch/FOLLOW_UP/NGC1851_metafiles/2022-08-02-06:50:21.meta" --candidate_csv '/b/PROCESSING/12_CANDYJAR/3HM_ACCELSEARCH_FULL/TOP_300/classifications/T1_CANDS.csv' --out_dir "/bscratch/CaReFuL/3HM_FOLLOW_UP/" --dm_half_range 0.2 --dm_tol 0.05    

if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description="Process beamforming metadata and candidate files.")
    parser.add_argument("--meta1", type=str, required=True, help="Path to the meta file for search observation.")
    parser.add_argument("--meta2", type=str, required=True, help="Path to the meta file for the follow up observation.")
    parser.add_argument("--candidate_csv", type=str, required=True, help="Path to the candidate CSV file.")
    parser.add_argument("--out_dir", type=str, required=True, help="Output directory for results.")
    parser.add_argument("--dm_half_range", type=float, default=0.2, help="Half range for DM values.")
    parser.add_argument("--dm_tol", type=float, default=0.05, help="Tolerance for DM values.")
    parser.add_argument("--run_name", type=str, default="3HM_FOLLOW_UP", help="Run name for the follow up.")
    parser.add_argument("--data_root", type=str, required=True, help="Root directory where data is kept. If this is Hercules archive, use sshfs to mount it to a local directory.")

    args = parser.parse_args()

    compact_meta = parse_meta_file(args.meta1)
    trapum_meta = parse_meta_file(args.meta2)
    candidate_df = pd.read_csv(args.candidate_csv)
    out_dir = args.out_dir + "/" + args.run_name
    data_root = args.data_root
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dm_half_range = args.dm_half_range
    dm_tol = args.dm_tol

    plt.clf()
    plt.figure(figsize=(15,15))
    ax = plt.gca()

    compact_beams = []
    trapum_beams = []
    colors = ['#3182bd', '#272727']
    line_styles = ['-', '--']
    alphas = [0.5, 0.3]


    for meta, color, line_style, alpha, beams, label in zip([compact_meta, trapum_meta], colors, line_styles, alphas, [compact_beams, trapum_beams], ["meta1", "meta2"]):

        boresight = meta["boresight"]
        target = boresight.split(",")[0].strip()
        pid = meta["project_name"]
        freq = int(float(meta["centre_frequency"]) /1e6)
        utc_start = meta["utc_start"].replace("/", "-").replace(" " ,"-")
        time = Time( meta["utc_start"].replace("/", "-").replace(" ", "T"), format='isot', scale='utc')

        bore_coords = SkyCoord(frame='icrs',
                            ra=boresight.split(",")[2].strip(),
                            dec=boresight.split(",")[3].strip(),
                            unit=(u.hourangle, u.deg))


        beamshape_json = json.loads(compact_meta["beamshape"])
        x = beamshape_json["x"]
        y =beamshape_json["y"]
        angle = beamshape_json["angle"]
        bore_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                bore_coords, bore_coords, time
            )
        bore_pixel_ra = bore_coords.ra.deg + bore_pixel_coordinates[0][0]
        bore_pixel_dec = bore_coords.dec.deg + bore_pixel_coordinates[0][1]




        for beam in meta["beams"]:
            if "ifbf" in beam:
                continue
            beam_values = meta["beams"][beam]
            if beam_values.split(",")[0].strip() == "unset":
                continue

            

            beam_coords = SkyCoord(frame='icrs',
                            ra=beam_values.split(",")[2].strip(),
                            dec=beam_values.split(",")[3].strip(),
                            unit=(u.hourangle, u.deg))
            pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                        beam_coords, bore_coords, time
                    )
            pixel_ra = bore_coords.ra.deg + pixel_coordinates[0][0]
            pixel_dec = bore_coords.dec.deg + pixel_coordinates[0][1]
            beams.append(Beam(beam, utc_start, pid, freq, target, beam_coords.ra.deg, beam_coords.dec.deg, pixel_ra, pixel_dec, x, y, angle))
            ellipse = Ellipse(xy=(pixel_ra, pixel_dec),
                                width=2*x, height=2*y,
                                angle=angle, facecolor='none', 
                                edgecolor=color, lw=1.5, alpha=alpha, linestyle=line_style)    
            ax.add_patch(ellipse)
            ax.plot(pixel_ra, pixel_dec, 'o', color=color, markersize=2, alpha=0.3)
            ax.annotate(beam.replace("3HM_", "").replace("cfbf00",""), (pixel_ra, pixel_dec), alpha=0.5, color=color, fontsize=8, ha='center', va='center')
    #plot unique legend
    #handles = [plt.Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor=color, markersize=10) for label, color in zip(labels, colors)]        
    extent = 0.02
    ax.set_title(f"meta1: {compact_meta['utc_start']} meta2: {trapum_meta['utc_start']}")
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    utc_dir = f"{out_dir}/{trapum_beams[0].utc}"
    if not os.path.exists(utc_dir):
        os.makedirs(utc_dir)
    plt.savefig(f"{utc_dir}/beams_full.pdf", bbox_inches='tight')
    for extent in [0.1, 0.075, 0.05, 0.025]:
        ax.set_xlim(bore_pixel_ra -extent, bore_pixel_ra + extent)
        ax.set_ylim(bore_pixel_dec - extent, bore_pixel_dec + extent)
        plt.savefig(f"{utc_dir}/beams_ex_{extent}.pdf", bbox_inches='tight')


    for compact_beam in compact_beams:
        overlapping_beams = [other_beam if discrete_overlap(compact_beam, other_beam) else None for other_beam in trapum_beams]
        #remove Nones if they exist
        overlapping_beams = [beam for beam in overlapping_beams if beam]
        compact_beam.add_overlapping_beams(overlapping_beams)
        print(compact_beam, ",".join([b.name for b in overlapping_beams]))


    #convert this into a map of rows for each beam_id value inside candidate_df so that I get a dataframe of candidates for each beam
    beam_candidate_map = {}
    for index, row in candidate_df.iterrows():
        beam_name = row["beam_name"]
        if beam_name not in beam_candidate_map:
            beam_candidate_map[beam_name] = []
        beam_candidate_map[beam_name].append(row)
        
    for beam_name in beam_candidate_map:
        compact_beam = [beam for beam in compact_beams if beam.name == beam_name][0]
        overlapping_beams = compact_beam.overlapping_beams
        candidate_df = pd.DataFrame(beam_candidate_map[beam_name])
        dm_list= []

        for dm in candidate_df['dm_opt'].values:
            dm_range = np.arange(dm - dm_half_range, dm + dm_half_range, dm_tol)
            dm_list.extend(dm_range)
        dm_list = sorted(set([f"{dm:.2f}" for dm in dm_list]))
        print(len(dm_list), dm_list)

        for overlapping_beam in overlapping_beams:
            beam_dir = f"{out_dir}/{overlapping_beam.utc}/{overlapping_beam.name}"
            if not os.path.exists(beam_dir):
                os.makedirs(beam_dir)
            dm_file = f"{out_dir}/{overlapping_beam.utc}/{overlapping_beam.name}/input.dmfile"
            with open(dm_file, 'w') as f:
                f.write("\n".join(dm_list))
            print(f"Written {dm_file}")

            cand_file = f"{out_dir}/{overlapping_beam.utc}/{overlapping_beam.name}/input.candfile"
            with open(cand_file, 'w') as cand_file_writer:
                cand_file_writer.write("#id DM accel F0 F1 F2 S/N\n")
                for i, row in candidate_df.iterrows():
                    cand_file_writer.write(f"{i} {row['dm_opt']} {row['acc_opt']} {row['f0_opt']} 0 0 {row['sn_fold']}\n")
            
            if data_root:
                data_dir_glob_str = f"{data_root}/*/*/{overlapping_beam.pid}/{overlapping_beam.target}/{overlapping_beam.utc}/{overlapping_beam.freq}"
                data_dirs = glob.glob(data_dir_glob_str)
                filterbanks = []
                for data_dir in data_dirs:
                    filterbanks.extend(glob.glob(f"{data_dir}/{beam}/*.fil"))
                
                filterbanks = sorted(set(filterbanks))
                print("Filterbanks: ", filterbanks)
                overlapping_beam.add_filterbanks(filterbanks)

                with open(f"{beam_dir}/input_fil.list", 'w') as filterbank_file:
                    for filterbank in filterbanks:
                        print(filterbank, data_root)
                        filterbank = filterbank.replace(data_root, "/tmp/UUID/")
                        filterbank_file.write(f"{filterbank}\n")
                print(f"Written {beam_dir}/input_fil.list")









