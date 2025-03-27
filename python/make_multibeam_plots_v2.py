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
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib.colors import Normalize
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from copy import copy

import argparse
params = {"figure.figsize": (12, 9),
          "font.size": 24,
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



candidate_root_dir="/bscratch/CAND_RECHECKS_2"
suffix="_fold-vals"
#cands_file_glob_pattern="*_no_*"

def get_logical_name(cands_file_glob_pattern):
    if "_no_" in cands_file_glob_pattern:
        return "no_reoptimisation"
    elif "__" in cands_file_glob_pattern:
        return "all_reoptimised"
    elif "nof1_nof0" in cands_file_glob_pattern:
        return "dm_reoptimised"
    else:
        return "Unknown"


cand_files = glob.glob(candidate_root_dir + "/*HM*.candfile")


df_map = {}
candidate_map = {}

for cand_file in cand_files:
    df = pd.read_csv(cand_file, sep="\s+")
    df_map[cand_file] = df
    with open(cand_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip()
            if line in candidate_map:
               candidate_map[line].append(cand_file)
            else:
                candidate_map[line] = [cand_file]

    

        
def get_non_comment_lines(file_path):
    with open(file_path) as f:
        return [line.strip() for line in f if not line.startswith("#")]

def get_index_of_line(file_path, key):
    try:
        return get_non_comment_lines(file_path).index(key)
    except ValueError:
        return -1
    
def get_line_of_index(file_path, index):
    lines = get_non_comment_lines(file_path)
    #remove lines with empty strings
    lines = [line for line in lines if line]
    return lines[index] if index < len(lines) else None


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

def pointInEllipse(x,y,xp,yp,d,D,angle):
    """
    tests if a point[xp,yp] is within
    boundaries defined by the ellipse
    of center[x,y], diameter d D, and tilted at angle
    """

    cosa=math.cos(angle)
    sina=math.sin(angle)
    dd=d/2*d/2
    DD=D/2*D/2

    a = math.pow(cosa*(xp-x)+sina*(yp-y),2)
    b = math.pow(sina*(xp-x)-cosa*(yp-y),2)
    ellipse=(a/dd)+(b/DD)
    if ellipse <= 1:
        return True
    else:
        return False

targets_file ="/b/PROCESSING/01_BEAMFORMED/J0514-4002A/2024-05-19-15:50:23/swdelays_J0514-4002A_1716133768_to_1716140969_f63918.targets"
out_dir="/bscratch/CAND_RECHECKS_2/output_pdfs"

targets = np.loadtxt(targets_file, delimiter=",", usecols=(0,1,2,3,4,5), dtype={
    'names': ('name', 'ra', 'dec', 'a', 'b', 'angle'),
    'formats': ('|S32', '|S32', '|S32', 'f8', 'f8', 'f8')
}, skiprows=1)

for cands_file_glob_pattern in ["*_no_*", "*__*", "*nof1_nof0*"]:

    candidate_beam_snrs = {}
    candiate_beam_duty_cycles = {}
    candidate_beam_plots = {}

    pdf_suffix = get_logical_name(cands_file_glob_pattern)

    print(f"Processing {pdf_suffix} candidates")

    for key, cand_file_list in candidate_map.items():
        print(f"Processing candidate {key}")
        beam_snr = {}
        beam_pngs = {}
        beam_duty_cycles = {}
        for cand_file in cand_file_list:
            beam = os.path.basename(cand_file).replace(".candfile", "").replace(suffix, "")
            candidate_idx = get_index_of_line(cand_file, key)
            if candidate_idx == -1:
                print("ERROR: Candidate not found in file")
                continue

            glob_pattern = os.path.join(candidate_root_dir, beam, beam + cands_file_glob_pattern + ".cands")
            cand_files_found = glob.glob(glob_pattern)
            if not cand_files_found:
                print(f"ERROR: No cands file for beam {beam}")
                continue
            candidate_cand_file = cand_files_found[0]


            candidate_line = get_line_of_index(candidate_cand_file, candidate_idx)
            if candidate_line is None or candidate_line == "":
                print("ERROR: Candidate number out of range")
                continue
            snr = candidate_line.split()[-2]
            width = candidate_line.split()[-1]
            f0 = candidate_line.split()[6]
            p0 = 1/float(f0)
            duty_cycle = 100 * float(width) / p0
            beam = "_".join(beam.split("_")[:-1])
            beam_snr[beam] = float(snr)
            beam_duty_cycles[beam] = duty_cycle
            glob_pattern = os.path.join(candidate_root_dir, os.path.basename(cand_file).replace(".candfile", "").replace(suffix, ""), f"*{beam}*_{cands_file_glob_pattern}{candidate_idx+1:04}.png")
            candidate_pngs = glob.glob(glob_pattern)
            if candidate_pngs:
                candidate_png = candidate_pngs[0]
                beam_pngs[beam] = candidate_png
            else:
                print(f"ERROR: No png found for candidate {key} and beam {beam} at {glob_pattern}")

        candidate_beam_snrs[key] = beam_snr
        candidate_beam_plots[key] = beam_pngs
        candiate_beam_duty_cycles[key] = beam_duty_cycles


    utc = '2024-05-19T15:50:23'
    ra = coord.Angle([x.decode() for x in targets['ra']], unit='hour')
    dec = coord.Angle([x.decode() for x in targets['dec']], unit='degree')
    names = targets['name']
    time = Time(utc, format='isot', scale='utc')

    boresight_ra_deg = ra[0].degree
    boresight_dec_deg = dec[0].degree

    bore_coords = SkyCoord(frame='icrs',
                        ra=boresight_ra_deg,
                        dec=boresight_dec_deg,
                        unit=(u.deg, u.deg))
    bore_pixel_coordinates = convert_equatorial_coordinate_to_pixel(
        bore_coords, bore_coords, time
    )
    bore_pixel_ra = boresight_ra_deg + bore_pixel_coordinates[0][0]
    bore_pixel_dec = boresight_dec_deg + bore_pixel_coordinates[0][1]

    import matplotlib.colors as mcolors
    color_dict = {
        '2HM': mcolors.CSS4_COLORS['slategray'],
        '3HM': mcolors.CSS4_COLORS['steelblue'],
        'T1':  mcolors.CSS4_COLORS['firebrick'],
        'T2':  mcolors.CSS4_COLORS['olive'],
    }

    names = np.unique([x.decode().split("_")[0] for x in names])

    #sort candidate_beam_snrs by candidate id
    candidate_beam_snrs = dict(sorted(candidate_beam_snrs.items(), key=lambda x: int(x[0].split()[0])))

    t1_cands_file="/homes/vkrishnan/dev/CaReFuL/python/input_sorted.csv"
    t1_cands = pd.read_csv(t1_cands_file, sep=",")
    count=0
    for key in candidate_beam_snrs.keys():

        beam_snrs_map = candidate_beam_snrs[key]
        beam_pngs_map = candidate_beam_plots[key]
        beam_duty_cycles_map = candiate_beam_duty_cycles[key]

        count+=1
        id = int(key.split()[0])
        P0 = 1000.0 / float(key.split()[3]) 
        DM = float(key.split()[1])

        pdf = PdfPages(os.path.join(out_dir, f'candidate_{id:03}_refolds_{pdf_suffix}.pdf'))

        #get idth row from t1_cands as index

        t1_cands_row = t1_cands.iloc[id]


        original_png_path = t1_cands_row['original_png_path']
        img = mpimg.imread(original_png_path)

        if not beam_snrs_map:
            continue

        beam_snrs = list(beam_snrs_map.values())
        beam_duty_cycles = list(beam_duty_cycles_map.values())
        beam_with_best_snr = max(beam_snrs_map, key=beam_snrs_map.get)
        # print(f"Best beam for candidate {id} is {beam_with_best_snr} with SNR {beam_snrs_map[beam_with_best_snr]}, "
        #       f"with range {min(beam_snrs)} to {max(beam_snrs)}")
        
        zoom_beam=beam_with_best_snr

        zoom_beam_target = [x for x in targets if x['name'].decode() == zoom_beam][0]
        zoom_beam_ra = coord.Angle( zoom_beam_target['ra'].decode(), unit='hour')
        zoom_beam_dec = coord.Angle( zoom_beam_target['dec'].decode(), unit='degree')
        zoom_beam_coords = SkyCoord(frame='icrs', ra=zoom_beam_ra.degree, dec=zoom_beam_dec.degree, unit=(u.deg, u.deg))
        zoom_pixel_coordinates = convert_equatorial_coordinate_to_pixel(zoom_beam_coords,bore_coords, time)
        zoom_pixel_beam_ra = boresight_ra_deg + zoom_pixel_coordinates[0][0]
        zoom_pixel_beam_dec = boresight_dec_deg + zoom_pixel_coordinates[0][1]

        #best 2HM beam snr 
        two_hm_beams = [x for x in beam_snrs_map.keys() if "2HM" in x]
        if two_hm_beams:
            best_2HM_beam = max(two_hm_beams, key=beam_snrs_map.get)
            best_2HM_snr = beam_snrs_map[best_2HM_beam]
            best_2HM_duty_cycle = beam_duty_cycles_map[best_2HM_beam]
        else:
            best_2HM_snr = 0
            best_2HM_duty_cycle = 0

        #best 3HM beam snr
        three_hm_beams = [x for x in beam_snrs_map.keys() if "3HM" in x]
        if three_hm_beams:
            best_3HM_beam = max(three_hm_beams, key=beam_snrs_map.get)
            best_3HM_snr = beam_snrs_map[best_3HM_beam]
            best_3HM_duty_cycle = beam_duty_cycles_map[best_3HM_beam]
        else:
            best_3HM_snr = 0
            best_3HM_duty_cycle = 0
        

        # Create figure and subplots
        plt.clf()
        # fig, axes = plt.subplots(figsize=(16, 9), nrows=2, ncols=2)
        # ax1, ax2, ax3, ax4 = axes.flatten()

        fig = plt.figure(figsize=(24,9))
        gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1], height_ratios=[1, 1])
        ax1 = plt.subplot(gs[0,0])
        ax2 = plt.subplot(gs[0,1])
        ax3 = plt.subplot(gs[1,0])
        ax4 = plt.subplot(gs[1,1])
        ax5 = plt.subplot(gs[:,2])
        axes = [ax1, ax2, ax3, ax4]

        axes_2hm = [ax1, ax2]
        axes_3hm = [ax3, ax4]

        axes_snr = [ax1, ax3]
        axes_duty_cycle = [ax2, ax4]

        ax1.set_title(f"Candidate {id}: Best 2HM S/N = {best_2HM_snr:.1f} ", fontsize=12)
        ax2.set_title(f"Candidate {id}: Best 2HM duty cycle = {best_2HM_duty_cycle:.1f}%  ", fontsize=12)
        ax3.set_title(f"Candidate {id}: Best 3HM S/N = {best_3HM_snr:.1f} ", fontsize=12)
        ax4.set_title(f"Candidate {id}: Best 3HM duty cycle = {best_3HM_duty_cycle:.1f}%  ", fontsize=12)

        # ax1.annotate("ax1", (0.5, 0.5), xycoords='axes fraction', va='center', ha='center', fontsize=50)
        # ax2.annotate("ax2", (0.5, 0.5), xycoords='axes fraction', va='center', ha='center', fontsize=50)
        # ax3.annotate("ax3", (0.5, 0.5), xycoords='axes fraction', va='center', ha='center', fontsize=50)
        # ax4.annotate("ax4", (0.5, 0.5), xycoords='axes fraction', va='center', ha='center', fontsize=50)
        # plt.savefig(f"test3.png")
        # break

        # Zoom extents
        extents_list = [0.015, 0.015,0.035, 0.035]
        for i, ax in enumerate(axes):
            extents = extents_list[i] 
            ax.set_ylabel('DEC (deg)', fontsize=14)
            ax.set_xlabel('RA (deg)', fontsize=14)
            ax.set_xlim(zoom_pixel_beam_ra - extents, zoom_pixel_beam_ra + extents)
            ax.set_ylim(zoom_pixel_beam_dec - extents, zoom_pixel_beam_dec + extents)


        # Setup color normalization and colormap
        snr_norm = Normalize(vmin=min(beam_snrs), vmax=max(beam_snrs))
        duty_cycle_norm = Normalize(vmin=min(beam_duty_cycles), vmax=max(beam_duty_cycles))
        snr_cmap = cm.coolwarm
        duty_cycle_cmap = cm.BrBG


        # Add ellipses for each beam
        for beam in beam_snrs_map.keys():
            snr = beam_snrs_map[beam]
            duty_cycle = beam_duty_cycles_map[beam]

            target = [x for x in targets if x['name'].decode() == beam][0]
            name = target['name'].decode()
            beam_ra = coord.Angle(target['ra'].decode(), unit='hour')
            beam_dec = coord.Angle(target['dec'].decode(), unit='degree')
            beam_width = 2 * target['a'].astype(float)
            beam_height = 2 * target['b'].astype(float)
            beam_angle = target['angle'].astype(float)
            beam_no = name.split("_")[1]

            # Convert equatorial beam coords to pixel coords
            beam_coords = SkyCoord(frame='icrs', ra=beam_ra.degree,
                                dec=beam_dec.degree, unit=(u.deg, u.deg))
            pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                beam_coords, bore_coords, time
            )
            pixel_beam_ra = boresight_ra_deg + pixel_coordinates[0][0]
            pixel_beam_dec = boresight_dec_deg + pixel_coordinates[0][1]

            colored_ellipse_snr = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec),
                            width=beam_width, height=beam_height,
                            angle=beam_angle, facecolor=snr_cmap(snr_norm(float(snr))),
                            edgecolor='black', lw=1.5, alpha=0.7)

            colored_ellipse_duty_cycle = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec),
                            width=beam_width, height=beam_height,
                            angle=beam_angle, facecolor=duty_cycle_cmap(duty_cycle_norm(float(duty_cycle))),
                            edgecolor='black', lw=1.5, alpha=0.7)
            
            hollow_ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec),
                                width=beam_width, height=beam_height,
                                angle=beam_angle, facecolor='none',
                                edgecolor='#272727', lw=1.5, alpha=0.4, linestyle='--')
            if '3HM' in name:
                for ax in axes_2hm:
                    ax.add_patch(copy(hollow_ellipse))
                    ax.annotate(f"{beam_no}",
                                    (pixel_beam_ra, pixel_beam_dec),
                                    fontsize=10, color="#272727", alpha=0.4)
                                    
                ax3.add_patch(colored_ellipse_snr)
                ax4.add_patch(colored_ellipse_duty_cycle)
                ax3.annotate(f"{beam_no}",
                                (pixel_beam_ra, pixel_beam_dec),
                                fontsize=10, color="black", alpha=1)
                ax4.annotate(f"{beam_no}",
                                (pixel_beam_ra, pixel_beam_dec),
                                fontsize=10, color="black", alpha=1)
            else:
            
                ax1.add_patch(colored_ellipse_snr)
                ax2.add_patch(colored_ellipse_duty_cycle)
                ax1.annotate(f"{beam_no}",
                                (pixel_beam_ra, pixel_beam_dec),
                                fontsize=10, color="black", alpha=1)
                ax2.annotate(f"{beam_no}",
                                (pixel_beam_ra, pixel_beam_dec),
                                fontsize=10, color="black", alpha=1)
                    


        # Adjust the layout to free space at the bottom (rect=...) or top
        # so the colorbar won't overlap the plots.
        plt.tight_layout(rect=[0, 0.08, 1, 1])  # Make room at the bottom

        for norm, axes, label, cmap in zip([snr_norm, duty_cycle_norm], [axes_snr, axes_duty_cycle], ["S/N", "Duty Cycle (%)"], [snr_cmap, duty_cycle_cmap]):

            # Create a single colorbar shared across both subplots
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  # Dummy array for the colorbar
            cbar = fig.colorbar(sm, ax=axes,
                                orientation="vertical",
                                fraction=0.05,   # Adjust fraction to control bar thickness
                                pad=0.08)       # Add padding so it's below the plots
            cbar.set_label(label, fontsize=12)
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.tick_params(labelsize=12)
            cbar.ax.yaxis.set_tick_params(width=1)
            cbar.ax.minorticks_on()



        #fig, ax = plt.subplots(figsize=(9,9))
        img = mpimg.imread(original_png_path)
        ax5.imshow(img, aspect='auto')
        plt.title(f"Original Candidate: {id}", fontsize=12)
        ax5.axis('off')
        pdf.savefig()
        plt.close(fig)

        sorted_beam_snrs = dict(sorted(beam_snrs_map.items(), key=lambda x: x[1], reverse=True))

        nfigs = len(sorted_beam_snrs.keys())//3 + 1
        for ifig in range(nfigs):
            fig, axes = plt.subplots(figsize=(24,9), ncols=3)   
            plotted = False 
            for col in range(3):
                ax = axes[col]
                ax.axis('off')
                idx = ifig*3 + col
                if idx >= len(sorted_beam_snrs.keys()):
                    continue
                beam = list(sorted_beam_snrs.keys())[idx]
                if beam not in beam_pngs_map:
                    print(f"ERROR: No png found for beam {beam}")
                    ax.annotate(f"No PNG found for beam {beam}", (0.5, 0.5), xycoords='axes fraction', va='center', ha='center', fontsize=12, color='red')
                    continue
                png = beam_pngs_map[beam]
                img = mpimg.imread(png)
                ax.imshow(img, aspect='auto')
                ax.set_title(f"{beam} S/N={sorted_beam_snrs[beam]:.1f}", fontsize=12)
                plotted = True
            if plotted:
                plt.tight_layout()
                pdf.savefig()
            plt.close(fig)

        plt.figure(figsize=(9,9))
        #pretty print t1_cands_row as a table in pdf
        ax = plt.gca()
        ax.axis('off')
        #save header and value in csv format 
        header = "Parameter,Value"
        values = [f"{k}:{v}" for k, v in t1_cands_row.items()]
        #introduce \n after every 150 characters
        formatted_values = []
        for value in values:
            #split value into lines of 80 characters
            chunks = [value[i:i+80] for i in range(0, len(value), 80)]
            formatted_values.extend(chunks)
            
        values = "\n".join(formatted_values)
        

        ax.text(0, 0, values, fontsize=12, ha='left', va='bottom')

        
        
        pdf.savefig()
        plt.close()
        pdf.close()



        
        #plt.savefig(f"candidate_{id}_beams.png")

        
        
