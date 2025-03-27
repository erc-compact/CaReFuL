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
candidate_root_dir="/bscratch/CAND_RECHECKS"
suffix="_fold-vals"
cands_file_glob_pattern="*_no_*"


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

targets = np.loadtxt(targets_file, delimiter=",", usecols=(0,1,2,3,4,5), dtype={
    'names': ('name', 'ra', 'dec', 'a', 'b', 'angle'),
    'formats': ('|S32', '|S32', '|S32', 'f8', 'f8', 'f8')
}, skiprows=1)


candidate_beam_snrs = {}
candidate_beam_plots = {}

for key, cand_file_list in candidate_map.items():
    beam_snr = {}
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
        if candidate_line is None:
            print("ERROR: Candidate number out of range")
            continue
        snr = candidate_line.split()[-1]
        beam = "_".join(beam.split("_")[:-1])
        beam_snr[beam] = float(snr)
        candidate_pngs = glob.glob(os.path.join(candidate_root_dir, os.path.basename(cand_file).replace(".candfile", "").replace(suffix, ""), f"*{beam}*_{cands_file_glob_pattern}{candidate_idx+1:04}.png"))
        if candidate_pngs:
            candidate_png = candidate_pngs[0]
            candidate_beam_plots[key] = candidate_png

    candidate_beam_snrs[key] = beam_snr
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
for key, beam_snrs_map in candidate_beam_snrs.items():
    count+=1
    if count > 1000:
        break
    id = int(key.split()[0])
    P0 = 1000.0 / float(key.split()[3]) 
    DM = float(key.split()[1])

    #get idth row from t1_cands as index

    t1_cands_row = t1_cands.iloc[id]

    print(t1_cands_row)

    original_png_path = t1_cands_row['original_png_path']
    img = mpimg.imread(original_png_path)

    if not beam_snrs_map:
        continue

    beam_snrs = list(beam_snrs_map.values())
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
    else:
        best_2HM_snr = 0

    #best 3HM beam snr
    three_hm_beams = [x for x in beam_snrs_map.keys() if "3HM" in x]
    if three_hm_beams:
        best_3HM_beam = max(three_hm_beams, key=beam_snrs_map.get)
        best_3HM_snr = beam_snrs_map[best_3HM_beam]
    else:
        best_3HM_snr = 0
    

    # Create figure and subplots
    plt.clf()
    fig = plt.figure(figsize=(16,9))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[:, 1])
    axes = [ax1, ax2, ax3]

    # Plot the original image
    ax3.imshow(img, aspect='auto')
    #ignore axis labels, ticks
    ax3.axis('off')

    #fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8,9))
                             #sharex=True, sharey=True)
    # axes[0].set_aspect(1)
    # axes[1].set_aspect(1)

    #axes[0].set_title(f"Candidate {id} P0={P0:.4f} ms DM={DM:.3f} pc/cc", fontsize=15)
    axes[0].set_title(f"Candidate {id} Best 2HM S/N = {best_2HM_snr:.1f} 3HM S/N = {best_3HM_snr:.1f}", fontsize=15)
    axes[0].set_ylabel('DEC (deg)', fontsize=14)
    axes[1].set_xlabel('RA (deg)', fontsize=14)
    axes[1].set_ylabel('DEC (deg)', fontsize=14)

    # Plot something minimal so each axes isn't empty
    axes[0].plot(bore_pixel_ra, bore_pixel_dec, '*', markersize=0.0)
    axes[1].plot(bore_pixel_ra, bore_pixel_dec, '*', markersize=0.0)

    # Zoom extents
    extents_list = [0.015, 0.035]
    for i, ax in enumerate(axes[:2]):
        extents = extents_list[i]
        ax.set_xlim(zoom_pixel_beam_ra - extents, zoom_pixel_beam_ra + extents)
        ax.set_ylim(zoom_pixel_beam_dec - extents, zoom_pixel_beam_dec + extents)

    # Setup color normalization and colormap
    norm = Normalize(vmin=min(beam_snrs), vmax=max(beam_snrs))
    cmap = cm.coolwarm

    # Add ellipses for each beam
    for beam, snr in beam_snrs_map.items():
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

        color = cmap(norm(float(snr)))
        colored_ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec),
                          width=beam_width, height=beam_height,
                          angle=beam_angle, facecolor=color,
                          edgecolor='black', lw=1.5, alpha=0.7)
        
        

        # Choose which of the 2 axes to draw on
        if "2HM" in name:
            ax_idx = 0
        else:
            hollow_ellipse = Ellipse(xy=(pixel_beam_ra, pixel_beam_dec),
                            width=beam_width, height=beam_height,
                            angle=beam_angle, facecolor='none',
                            edgecolor='#272727', lw=1.5, alpha=0.4, linestyle='--')
            axes[0].add_patch(hollow_ellipse)
            axes[ax_idx].annotate(f"{beam_no}",
                              (pixel_beam_ra, pixel_beam_dec),
                              fontsize=10, color="#272727", alpha=0.4)

            ax_idx = 1

        
        axes[ax_idx].add_patch(colored_ellipse)
        axes[ax_idx].annotate(f"{beam_no}",
                              (pixel_beam_ra, pixel_beam_dec),
                              fontsize=10, color="black", alpha=1)

    # Adjust the layout to free space at the bottom (rect=...) or top
    # so the colorbar won't overlap the plots.
    plt.tight_layout(rect=[0, 0.08, 1, 1])  # Make room at the bottom

    # Create a single colorbar shared across both subplots
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Dummy array for the colorbar
    cbar = fig.colorbar(sm, ax=axes[:2],
                        orientation="vertical",
                        fraction=0.05,   # Adjust fraction to control bar thickness
                        pad=0.08)       # Add padding so it's below the plots
    cbar.set_label("SNR", size=14)
    cbar.ax.tick_params(labelsize=12)
    #add minor ticks to colorbar
    cbar.ax.yaxis.set_tick_params(width=1)
    cbar.ax.minorticks_on()
    
    plt.savefig(f"candidate_{id}_beams.png")
    plt.close(fig)
    
    
