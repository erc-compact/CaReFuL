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
    def __init__(self, name, utc, ra, dec,pixel_ra, pixel_dec, x, y, angle):
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
    
    def add_overlapping_beams(self, beams):
        self.overlapping_beams.extend(beams)
    
    def __str__(self):
        return f"Beam {self.name} "
    
    def __repr__(self):
        return f"Beam {self.name} "
    
    def __eq__(self, other):
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)
    
    

compact_meta = parse_meta_file("/b/PROCESSING/01_BEAMFORMED/J0514-4002A/2024-05-19-15:50:23/3HM.meta")
trapum_meta = parse_meta_file("/bscratch/FOLLOW_UP/NGC1851_metafiles/2022-08-02-06:50:21.meta")
candidate_df = pd.read_csv('/b/PROCESSING/12_CANDYJAR/3HM_ACCELSEARCH_FULL/TOP_300/classifications/T1_CANDS.csv')
out_dir = "/bscratch/CaReFuL/3HM_FOLLOW_UP/"
dm_half_range = 0.2
dm_tol = 0.05

compact_beams = []
trapum_beams = []

for meta, color, beams in zip([compact_meta, trapum_meta], ['red', 'blue'], [compact_beams, trapum_beams]):
    boresight = meta["boresight"]
    utc_start = meta["utc_start"].replace("/", "-").replace(" ", "T")
    time = Time(utc_start, format='isot', scale='utc')

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
        beam_coords = SkyCoord(frame='icrs',
                        ra=beam_values.split(",")[2].strip(),
                        dec=beam_values.split(",")[3].strip(),
                        unit=(u.hourangle, u.deg))
        pixel_coordinates = convert_equatorial_coordinate_to_pixel(
                    beam_coords, bore_coords, time
                )
        pixel_ra = bore_coords.ra.deg + pixel_coordinates[0][0]
        pixel_dec = bore_coords.dec.deg + pixel_coordinates[0][1]
        beams.append(Beam(beam, utc_start, beam_coords.ra.deg, beam_coords.dec.deg, pixel_ra, pixel_dec, x, y, angle))


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
        dm_file = f"{out_dir}/{overlapping_beam.utc}/{overlapping_beam.name}/dm.file"
        with open(dm_file, 'w') as f:
            f.write("\n".join(dm_list))
        print(f"Written {dm_file}")

        cand_file = f"{out_dir}/{overlapping_beam.utc}/{overlapping_beam.name}/input.candfile"
        with open(cand_file, 'w') as cand_file_writer:
            cand_file_writer.write("#id DM accel F0 F1 F2 S/N\n")
            for i, row in candidate_df.iterrows():
                cand_file_writer.write(f"{i} {row['dm_opt']} {row['acc_opt']} {row['f0_opt']} 0 0 {row['sn_fold']}\n")

    
        







