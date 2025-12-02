"""
Run a synthpop model and turn the output into simulated images.

Instructions:

Edit the field output_dir_root in the script below to hard-wire in the
path to your synthpop output directory, i.e. the directory you have
initially set up with

> python -m synthpop.migrate_interactive_part path_to_directory

You may need to also edit the sys.path.append line so that this script
can find your synthpop installation directory.

Prepare a synthpop json configuration file similar to

{   "SEED":{"random_seed":null},

    "MANDATORY":{
        "#comment1": "directory and base for the output files",
        "model_name":"Huston2025",
        "#comment2": "directory containing population json files",
        "name_for_output":"Huston2025"
    },

    "SIGHTLINES":
        {
            "l_set": [0.5], "l_set_type":"list",
            "b_set":[1.5], "b_set_type":"list",
            "solid_angle": 2.7e-2, "solid_angle_unit": "deg^2"
        },

    "EXTINCTION_MAP":
        {
        "extinction_map_kwargs": {"name":"Surot", "project_3d":true, "dist_2d":8.15},
        "extinction_law_kwargs": [{"name":"SODC", "R_V":2.5}]
        },

    "POPULATION_GENERATION":{
        "skip_lowmass_stars": false
    },

    "PHOTOMETRIC_OUTPUTS":{
        "maglim":["W146", 99, "keep"],
        "chosen_bands": ["R062","Z087","Y106","J129","W146","H158","F184", "Bessell_U", "Bessell_B", "Bessell_V", "Bessell_R", "Bessell_I", "VISTA_J", "VISTA_H", "VISTA_Ks"]
    },

    "OUTPUT":{
        "post_processing_kwargs": [{"name":"ProcessDarkCompactObjects", "remove":false},
                {"name":"ConvertMistMags", "conversions":{"AB": ["R062", "Z087", "Y106", "J129", "W146", "H158", "F184"]}},
                {"name":"RenameColumns",
                    "old_names":["log_L", "log_Teff", "log_g", "[Fe/H]","log_R"],
                    "new_names":["logL", "logTeff", "logg" ,"Fe/H_evolved","log_radius"]}],

        "output_location":"outputfiles/lens",
        "output_filename_pattern": "{name_for_output}_l{l_deg:.3f}_b{b_deg:.3f}",

        "overwrite": true
    },

    "IMAGES": {
        "delta_t_minutes": [2, 15]
        },

    "MICROLENSING": {
        "microlensing_event_parameters": {
            "mag_bin_edges": [20, 21, 22, 23, 24, 25, 26]},
            "events_u0_per_mag_bin": [[200, 0.1], [200, 0.01], [200, 0.001]],
            "t_E_mins": 400,
            "t0_mins": 2160
        }
}

This is more-or-less the same as a standard synthpop config file, except only
the first l_set and b_set elements will be considered, and you need to provide
a delta_t_minutes specification of the number of images and their cadence.
One output image is made for each epoch.

The final field, "microlensing_event_parameters", if provided, is used to configure
PSPL microlensing events to be injected into random stars.

If your config file is config.json, Run this script with python make_synthpop_image.py config.json

"""

__author__ = "Michael Albrow"

# Point this to your synthpop output directory
# For kerr
output_dir_root = '/home/users/mda45/local/data/synthpop/outputfiles/'
# For rch
#output_dir_root = '/home/mda45/synthpop/synthpop_data/outputfiles/'




# Modules we will need
import sys
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
from copy import deepcopy
from functools import partial
from contextlib import redirect_stdout
import numpy as np
import pandas as pd
import argparse
from multiprocessing import Pool
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactic, FK5
import json
import synthpop
from galsim import UniformDeviate
from romanisim import wcs, persistence, parameters
from romanisim import ris_make_utils as ris
from asdf_to_fits import asdf_to_fits

MAX_PARALLEL_PROCESSES = int(os.cpu_count() / 2)

# These offsets shift the synthpop field to the approximate centre of SCA 1.
COORDINATE_OFFSET_RA_DEG = 0.0655
COORDINATE_OFFSET_DEC_DEG = 0.0459

ROMAN_PIXEL_ARCSEC_PER_PIXEL = 0.11

MINUTES_PER_DAY = 60 * 24
DAYS_PER_YEAR = 365.25
MINUTES_PER_YEAR = MINUTES_PER_DAY * DAYS_PER_YEAR

DEFAULT_BANDPASS = "F146"
DEFAULT_SCA = 1
DEFAULT_DATE = "2000-03-30T00:00:00"
DEFAULT_MA_TABLE_NUMBER = 4


def synthpop_to_romanisim(t: table.Table, delta_t_years: float,
                          out_file: (type(None), str) = None) -> table.Table:
    """
    Convert synthpop table, t, into romanisim required format.
    Write it to out_file if provided.
    The input table must at least have columns labelled 'l', 'b', 'mul', 'mub', and 'W146'.
    A 2D random dither is applied, as drawn from a normal distribution with the indicated amplitude.
    Proper motions (mul, mub) (in mas) are added to the coordinates using the supplied delta_t_years.
    Returns the catalogue as an astropy table.
    """

    c = SkyCoord(l=t['l'] * u.degree,
                 b=t['b'] * u.degree,
                 distance=t['Dist'] * u.kpc,
                 pm_l_cosb=t['mul'] * u.mas / u.year,
                 pm_b=t['mub'] * u.mas / u.year,
                 frame=Galactic)

    c.apply_space_motion(dt=delta_t_years * u.year)
    c_cel = c.fk5

    new_cat = table.Table([c_cel.ra, c_cel.dec], names=('ra', 'dec'), dtype=[np.float64, np.float64])

    new_cat['ra'].unit = None
    new_cat['dec'].unit = None
    new_cat['type'] = 'PSF'
    new_cat['n'] = -1.0
    new_cat['half_light_radius'] = 0.0
    new_cat['pa'] = 0.0
    new_cat['ba'] = 1.0
    new_cat['F146'] = 10.0 ** (-0.4 * t['W146'])

    if out_file is not None:
        print(f'Writing {out_file}')
        new_cat.write(out_file, format='ascii.ecsv', comment='#', delimiter=' ', overwrite=True)

    return new_cat


def configure_microlensing_event_parameters(parameters: dict, stars_table: table.Table,
                                            out_file='ulens_stars') -> dict:
    """Select stars for microlensing and add to parameters dict."""

    microlens_parameters = deepcopy(parameters)
    microlens_parameters['u0_stars'] = []

    cat = synthpop_to_romanisim(stars_table, 0.0)

    try:
        ulens_stars = np.loadtxt(out_file)
        for line in ulens_stars:
            microlens_parameters['u0_stars'].append([line[4], int(line[0])])

    except FileNotFoundError:

        print(f"ulens file {out_file} not found. Creating ...")

        with open(out_file, 'w') as f:

            for m1, m2 in zip(parameters['mag_bin_edges'][:-1], parameters['mag_bin_edges'][1:]):

                p = np.where((m1 < stars_table['W146']) & (stars_table['W146'] <= m2))[0]

                print(f'{len(p)} stars with mags between {m1} and {m2}')

                print('events_u0_per_mag_bin', parameters['events_u0_per_mag_bin'])
                for n_events, u0 in parameters['events_u0_per_mag_bin']:
                    size_select = min(n_events, len(p))
                    print(f'Selecting {size_select} events')

                    if size_select > 0:
                        p_select = np.random.choice(p, size=n_events, replace=False)
                        microlens_parameters['u0_stars'].append([u0, p_select])
                        for psel in p_select:
                            ra = cat['ra'][psel]
                            dec = cat['dec'][psel]
                            mag = stars_table['W146'][psel]
                            print(f'{psel}: {ra}, {dec}, {mag} {u0}')
                            f.write(f'{psel} {ra} {dec} {mag} {u0}\n')
                    else:
                        microlens_parameters["u0_stars"].append([u0, []])

    return microlens_parameters


def insert_microlensing_events(cat: table.table, t: float, event_params: dict) -> table.table:
    """Insert microlensing events into selected stars in the catalogue."""

    t_mins = t * MINUTES_PER_YEAR
    tau2 = ((t_mins - event_params['t0_mins']) / event_params['t_E_mins'])**2

    for u0, stars in event_params['u0_stars']:
        u = np.sqrt(u0**2 + tau2)
        for star in np.atleast_1d(stars):
            cat[star]['F146'] = cat[star]['F146'] * (u**2 + 2.0) / (u*np.sqrt(u**2 + 4.0))

    return cat


def simulate_args(ra: float, dec: float, filename: str) -> argparse.Namespace:
    """Build a minimal parser namespace needed by ris.simulate_image_file.

    This is mostly a copy of code from the romanisim-make-image script.
    """

    parser = argparse.ArgumentParser(
        description='Make a demo image.',
        epilog='EXAMPLE: %(prog)s output_image.asdf',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('--filename', type=str, help='output image (asdf)')
    parser.add_argument('--bandpass', type=str, help='bandpass to simulate',
                        default='F087')
    parser.add_argument('--boresight', action='store_true', default=False,
                        help=('radec specifies location of boresight, not '
                              'center of WFI.'))
    parser.add_argument('--catalog', type=str, help='input catalog (ecsv)',
                        default=None)
    parser.add_argument('--config', type=str, help='input parameter override file (yaml)',
                        default=None)
    parser.add_argument('--date', type=str, default=None,
                        help='UTC Date and Time of observation to simulate in ISOT format.')
    parser.add_argument('--level', type=int, default=2,
                        help='1 or 2, for L1 or L2 output')
    parser.add_argument('--ma_table_number', type=int, default=DEFAULT_MA_TABLE_NUMBER)
    parser.add_argument('--nobj', type=int, default=1000)
    parser.add_argument('--previous', default=None, type=str,
                        help=('previous simulated file in chronological order '
                              'used for persistence modeling.'))
    parser.add_argument('--radec', type=float, nargs=2,
                        help='ra and dec (deg)', default=None)
    parser.add_argument('--rng_seed', type=int, default=None)
    parser.add_argument('--roll', type=float, default=0,
                        help='Position angle (North towards YIdl) measured at the V2Ref/V3Ref of the aperture used.')
    parser.add_argument('--sca', type=int, default=7, help='SCA to simulate')
    parser.add_argument('--usecrds', action='store_true',
                        help='Use CRDS for distortion map')
    parser.add_argument('--stpsf', action='store_true',
                        help='Use stpsf for PSF')
    parser.add_argument('--truncate', type=int, default=None, help=(
        'If set, truncate the MA table at given number of resultants.'))
    parser.add_argument('--pretend-spectral', type=str, default=None, help=(
        'Pretend the image is spectral.  exposure.type and instrument.element '
        'are updated to be grism / prism.'))
    parser.add_argument('--drop-extra-dq', default=False, action='store_true',
                        help=('Do not store the optional simulated dq array.'))
    parser.add_argument('--scale-factor', type=float, default=-1.,
                        help=(
                            'Velocity aberration-induced scale factor. If negative, use given time to calculated based on orbit ephemeris.'))

    args = parser.parse_args(['--radec', f'{ra}', f'{dec}',
                              '--stpsf',
                              '--filename', filename,
                              '--date', DEFAULT_DATE,
                              '--sca', f"{DEFAULT_SCA}",
                              '--bandpass', DEFAULT_BANDPASS])

    return args


def dither(dither_pattern: str = "random", i: int = 0) -> (float, float):
    "Return an (x, y) subpixel dither based on the chosen pattern."

    if dither_pattern == "random":
        np.random.seed()
        return (np.random.rand() - 0.5, np.random.rand() - 0.5)

    elif dither_pattern == "Anderson_8x8":
        dx = [0.0, 4.5, 0.0, 4.5, 2.2, 6.7, 2.7, 6.7, 0.0, 4.5, 0.0, 4.5, 2.2, 6.7, 2.2, 6.7, \
              1.1, 5.6, 1.1, 5.6, 3.3, 7.8, 3.3, 7.8, 1.1, 5.6, 1.1, 5.6, 3.3, 7.8, 3.3, 7.8, \
              0.0, 4.5, 0.0, 4.5, 2.2, 6.7, 2.7, 6.7, 0.0, 4.5, 0.0, 4.5, 2.2, 6.7, 2.2, 6.7, \
              1.1, 5.6, 1.1, 5.6, 3.3, 7.8, 3.3, 7.8, 1.1, 5.6, 1.1, 5.6, 3.3, 7.8, 3.3, 7.8]
        dy = [0.0, 0.0, 4.5, 4.5, 0.0, 0.0, 4.5, 4.5, 2.2, 2.2, 6.7, 6.7, 2.2, 2.2, 6.7, 6.7, \
              0.0, 0.0, 4.5, 4.5, 0.0, 0.0, 4.5, 4.5, 2.2, 2.2, 6.7, 6.7, 2.2, 2.2, 6.7, 6.7, \
              1.1, 1.1, 5.6, 5.6, 1.1, 1.1, 5.6, 5.6, 3.3, 3.3, 7.8, 7.8, 3.3, 3.3, 7.8, 7.8, \
              1.1, 1.1, 5.6, 5.6, 1.1, 1.1, 5.6, 5.6, 3.3, 3.3, 7.8, 7.8, 3.3, 3.3, 7.8, 7.8]
        dxi = dx[i] % 64
        dyi = dy[i] % 64
        return dxi, dyi

    raise ValueError(f'dither pattern {dither_pattern} not recognized')


def make_image(i: int, delta_t: float, file_root: str, star_table: table.Table, ra: float, dec: float,
               dither_pattern: str = "random", microlensing_event_parameters: dict = None) -> None:
    """Make a single image with romanisim."""

    with open(f'{file_root}.log', 'w', buffering=1) as f:
        with redirect_stdout(f):

            cat = synthpop_to_romanisim(star_table, delta_t)

            if microlensing_event_parameters is not None:
                cat = insert_microlensing_events(cat, delta_t, microlensing_event_parameters)

            # Random dither
            roman_pixel_scale = u.pixel_scale(ROMAN_PIXEL_ARCSEC_PER_PIXEL * u.arcsec / u.pixel)
            dx_pixels, dy_pixels = dither(dither_pattern)
            d_ra = (dx_pixels * u.pixel).to(u.degree, roman_pixel_scale) / np.cos(dec * u.degree)
            d_dec = (dy_pixels * u.pixel).to(u.degree, roman_pixel_scale)
            print(f"dither: ({dx_pixels}, {dy_pixels}) pixels,   ({d_ra.value}, {d_dec.value}) degrees")

            # Create persistence object - not needed?
            persist = persistence.Persistence()

            args = simulate_args(ra, dec, f'{file_root}.asdf')

            metadata = ris.set_metadata(
                date=args.date, bandpass=args.bandpass,
                sca=args.sca, ma_table_number=args.ma_table_number,
                truncate=args.truncate)

            coord = SkyCoord(ra=args.radec[0] * u.deg + d_ra, dec=args.radec[1] * u.deg + d_dec, frame='icrs')
            wcs.fill_in_parameters(metadata, coord, boresight=args.boresight, pa_aper=args.roll)

            # Simulate image and write to file in asdf format
            rng = UniformDeviate(None)
            ris.simulate_image_file(args, metadata, cat, rng, persist)

            # Also write the image in FITS format
            asdf_to_fits(f'{file_root}.asdf', f'{file_root}.fits')

            sys.stdout.flush()


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Usage: python make_synthpop_image.py config.json")
        sys.exit(1)

    config_file = sys.argv[1]

    try:
        with open(config_file) as file:
            config_data = json.load(file)
    except FileNotFoundError:
        print(f"Config file {config_file} not found.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Invalid JSON in config file {config_file}.")
        sys.exit(1)

    # Centre-of-field coordinates for romanisim
    gal_l = config_data["SIGHTLINES"]["l_set"][0]
    gal_b = config_data["SIGHTLINES"]["b_set"][0]
    c = SkyCoord(l=gal_l * u.degree, b=gal_b * u.degree, frame=Galactic)
    c_cel = c.fk5
    ra_rom = c_cel.ra + COORDINATE_OFFSET_RA_DEG * u.degree
    dec_rom = c_cel.dec + COORDINATE_OFFSET_DEC_DEG * u.degree

    # Run synthpop if catalogue doesn't yet exist
    synthpop_cat_file = \
        f'{config_data["OUTPUT"]["output_location"]}/{config_data["MANDATORY"]["model_name"]}_l{gal_l:.3f}_b{gal_b:.3f}.csv'
    csv_used_columns = ['W146', 'l', 'b', 'Dist', 'mul', 'mub']
    try:
        df = pd.read_csv(synthpop_cat_file, usecols=csv_used_columns)
    except FileNotFoundError:
        model = synthpop.SynthPop(config_file, overwrite=True)
        model.init_populations()
        model.process_all()
        print('synthpop output_location:', model.parms.output_location)
        df = pd.read_csv(synthpop_cat_file, usecols=csv_used_columns)

    # Set offset time epochs
    n_images, cadence_minutes = config_data["IMAGES"]["delta_t_minutes"]
    t_minutes = np.arange(n_images) * cadence_minutes
    t_years = t_minutes / MINUTES_PER_YEAR

    # Read synthpop output catalogue
    synthpop_table = table.Table.from_pandas(df)
    synthpop_table = synthpop_table[~synthpop_table['W146'].mask]

    # Configure microlensing events to insert
    ulens_parameters = configure_microlensing_event_parameters(
        config_data["MICROLENSING"]["microlensing_event_parameters"],
        synthpop_table)
    try:
        ulens_parameters = configure_microlensing_event_parameters(config_data["MICROLENSING"]["microlensing_event_parameters"],
                                                                   synthpop_table)
    except KeyError:
        ulens_parameters = None

    print("ulens_parameters:", ulens_parameters)

    if "dither_pattern" not in config_data["IMAGES"].keys():
        dither_pattern = "random"
    else:
        dither_pattern = config_data["IMAGES"]["dither_pattern"]

    if "output_location" not in config_data["IMAGES"].keys():
        images_output_location = "."
    else:
        images_output_location = config_data["IMAGES"]["output_location"]

    # Make images
    file_roots = [f'{images_output_location}/{config_data["MANDATORY"]["name_for_output"]}_t{i:04d}' for i in range(n_images)]

    if MAX_PARALLEL_PROCESSES > 1 and n_images > 1:

        n_processes = min(MAX_PARALLEL_PROCESSES, n_images)

        with Pool(n_processes) as pool:
            pool.starmap(partial(make_image, star_table=synthpop_table, ra=ra_rom.value, dec=dec_rom.value,
                                 dither_pattern=dither_pattern, microlensing_event_parameters=ulens_parameters),
                         zip(range(len(t_years)), t_years, file_roots))

    else:

        for i, (t, file_root) in enumerate(zip(t_years, file_roots)):
            make_image(i, t, file_root, synthpop_table, ra_rom.value, dec_rom.value,
                       dither_pattern=dither_pattern, microlensing_event_parameters=ulens_parameters)

