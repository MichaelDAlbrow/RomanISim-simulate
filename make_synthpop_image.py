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

{
    "model_base_name" : "my_generated_model",
    "l_set" : [0],
    "b_set" : [-1],
    "solid_angle" : 3e-3,
    "model_name" : "besancon_Robin2003",
    "evolution_class" : {"name" : "MIST", "interpolator" : "CharonInterpolator"},
    "extinction_map_kwargs" : {"name":"Marshall"},
    "extinction_law_kwargs" : {"name":"ODonnellCardelli"},
    "delta_t_minutes" :  [100, 15],
    "microlensing_event_parameters" : {
        "mag_bin_edges" : [20, 21, 22, 23, 24, 25, 26],
        "events_u0_per_mag_bin" : [[20, 0.1], [20, 0.01]],
        "t_E_mins" : 288,
        "t0_mins" : 2160
  }
}

This is more-or-less the same as a standard synthpop config file, except only the
first l_set and b_set elements will be considered, and you need to provide a
delta_t_minutes specification of the number of images and their cadence.

The final field, "microlensing_event_parameters", if provided, is used to configure
PSPL microlensing events to be injected into random stars.

One output image is made for each epoch.

If your config file is config.json, Run this script with
python make_synthpop_image.py config.json

"""

__author__ = "Michael Albrow"

# Point this to your synthpop output directory
output_dir_root = '/scratch/mda45/synthpop-data/outputfiles/'

# Modules we will need
import sys
import os
from copy import deepcopy
from functools import partial
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
from romanisim import log, wcs, persistence, parameters
from romanisim import ris_make_utils as ris
from asdf_to_fits import asdf_to_fits

# To profile the romanISIM code
# import cProfile

#  Uncomment and edit if necessary
sys.path.append('/Users/mda45/Packages/synthpop-main')

# For parallel processing. Set this to 1 if you don't want parallel processing.
# max_parallel_processes = 1
max_parallel_processes = np.min([200, os.cpu_count()])

# These offsets shift the synthpop field to the approximate centre of SCA 1.
coordinate_offset = {'ra': 0.0655 * u.degree, 'dec': 0.0459 * u.degree}


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

    cat = synthpop_to_romanisim(stars_table, 0.0)

    with open(out_file, 'w') as f:

        for m1, m2 in zip(parameters['mag_bin_edges'][:-1], parameters['mag_bin_edges'][1:]):

            p = np.where((m1 < stars_table['W146']) & (stars_table['W146'] <= m2))[0]

            microlens_parameters['u0_stars'] = []
            for n_events, u0 in parameters['events_u0_per_mag_bin']:
                p_select = np.random.choice(p, size=n_events, replace=False)
                microlens_parameters['u0_stars'].append([u0, p_select])

            for psel in p_select:
                ra = cat['ra'][psel]
                dec = cat['dec'][psel]
                mag = stars_table['W146'][psel]
                f.write(f'{psel} {ra} {dec} {mag} {u0}\n')

    return microlens_parameters


def insert_microlensing_events(cat: table.table, t: float, event_params: dict) -> table.table:
    """Insert microlensing events into selected stars in the catalogue."""

    tau2 = ((t - event_params['t0_mins']) / event_params['t_E_mins'])**2

    for u0, stars in event_params['u0_stars']:
        u = np.sqrt(u0**2 + tau2)
        for star in stars:
            cat[star]['F146'] = cat[star]['F146'] * (u**2 + 2.0) / (u*np.sqrt(u**2 + 4.0))

    return cat


def simulate_args(ra: float, dec: float, filename: str) -> argparse.Namespace:
    """Build a minimal parser namespace needed by ris.simulate_image_file.

    This is mostly a copy of code from the romanisim-make-image script.
    """

    parser = argparse.ArgumentParser(
        description='Make a demo image.',
        epilog='EXAMPLE: %(prog)s output_image.asdf',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('--ma_table_number', type=int, default=1)
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
    parser.add_argument('--webbpsf', action='store_true',
                        help='Use webbpsf for PSF')
    parser.add_argument('--truncate', type=int, default=None, help=(
        'If set, truncate the MA table at given number of resultants.'))

    args = parser.parse_args(['--radec', f'{ra}', f'{dec}',
                              '--webbpsf',
                              '--filename', filename,
                              '--date', '2000-03-30T00:00:00',
                              '--sca', '1',
                              '--bandpass', 'F146'])

    return args


def make_image(delta_t: float, file_root: str, star_table: table.Table, ra: float, dec: float,
               dither_rms_pixels: float = 0.0, microlensing_event_parameters: (type(None) | dict) = None) -> None:
    """Make a single image with romanisim."""

    cat = synthpop_to_romanisim(star_table, delta_t)

    if microlensing_event_parameters is not None:
        cat = insert_microlensing_events(cat, delta_t, microlensing_event_parameters)

    rng = UniformDeviate(None)

    # Random dither
    roman_pixel_scale = u.pixel_scale(0.11 * u.arcsec / u.pixel)
    np.random.seed()
    dx_pixels, dy_pixels = dither_rms_pixels * np.random.randn(2)
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
    ris.simulate_image_file(args, metadata, cat, rng, persist)

    # Also write the image in FITS format
    asdf_to_fits(f'{file_root}.asdf', f'{file_root}.fits')


if __name__ == '__main__':

    config_file = sys.argv[1]

    with open(config_file) as file:
        config_data = json.load(file)

    # Run synthpop
    model = synthpop.SynthPop(config_file, overwrite=True)
    model.init_populations()
    model.process_all()
    print('synthpop output_location:', model.parms.output_location)

    # Set offset time epochs
    n_images, cadence_minutes = config_data["delta_t_minutes"]
    t_minutes = np.arange(n_images) * cadence_minutes
    t_years = t_minutes / (60 * 24 * 365.25)

    # Centre-of-field coordinates for romanisim
    gal_l = config_data["l_set"][0]
    gal_b = config_data["b_set"][0]
    c = SkyCoord(l=gal_l * u.degree, b=gal_b * u.degree, frame=Galactic)
    c_cel = c.fk5
    ra_rom = c_cel.ra + coordinate_offset['ra']
    dec_rom = c_cel.dec + coordinate_offset['dec']

    # Read synthpop output catalogue
    synthpop_cat_file = \
        f'{output_dir_root}/{config_data["name_for_output"]}/{config_data["model_name"]}_l{gal_l:.3f}_b{gal_b:.3f}.csv'
    df = pd.read_csv(synthpop_cat_file)
    synthpop_table = table.Table.from_pandas(df)
    synthpop_table = synthpop_table[~synthpop_table['W146'].mask]

    # Configure microlensing events to insert
    if "microlensing_event_parameters" in config_data:
        ulens_parameters = configure_microlensing_event_parameters(config_data["microlensing_event_parameters"],
                                                                   synthpop_table)
    else:
        ulens_parameters = None

    # Make images
    if max_parallel_processes > 1 and n_images > 1:

        with Pool(max_parallel_processes) as pool:
            pool.starmap(partial(make_image, star_table=synthpop_table, ra=ra_rom.value, dec=dec_rom.value,
                                 dither_rms_pixels=1.0, microlensing_event_parameters=ulens_parameters),
                         zip(t_years, [f'{config_data["name_for_output"]}_t{i:04d}' for i in range(n_images)]))

    else:

        # Uncomment for code profiling
        # pr = cProfile.Profile()
        # pr.enable()

        for i, t in enumerate(t_years):
            make_image(t, f'{config_data["name_for_output"]}_t{i:04d}', synthpop_table, ra_rom.value, dec_rom.value,
                       dither_rms_pixels=1.0, microlensing_event_parameters=ulens_parameters)

        # Uncomment for code profiling
        # pr.disable()
        # pr.print_stats(sort='cumtime')
