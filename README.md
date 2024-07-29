# RomanISim-simulate

Simulate images for the Nancy Grace Roman Space Telescope Time-Domain Survey.

This code uses synthpop and RomanISIM/WebbPSF to generate time-series simulated
images from Roman. 

A star catalog is first generated using synthpop. Then, a series of images (for detector SCA1)
are generated at a given cadence, applying the stellar kinematics from the catalog, and 
random dithering. 

Optionally, PSPL microlensing events can be inserted onto random source stars.


## Requires

python 3.9

## Required external packages

Install and configure https://github.com/synthpop-galaxy/synthpop.

Install and configure https://github.com/spacetelescope/romanisim (slightly  non-trivial) 
and https://github.com/spacetelescope/webbpsf.


## Run example

python -u make_synthpop_image.py synthpop_config18.json >& test18.log &

## More-detailed instructions

Edit the field output_dir_root in the make_synthpop_image.py to hard-wire in the
path to your synthpop output directory, i.e. the directory you have
initially set up with

> python -m synthpop.migrate_interactive_part path_to_directory

You may need to also edit the sys.path.append line so that the script
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
    "delta_t_minutes" :  [192, 15],
    "microlensing_event_parameters" : {
        "mag_bin_edges" : [20, 21, 22, 23, 24, 25, 26],
        "events_u0_per_mag_bin" : [[40, 0.1], [40, 0.01], [40, 0.001]],
        "t_E_mins" : 288,
        "t0_mins" : 2160
  }
}

This is more-or-less the same as a standard synthpop config file, except only the
first l_set and b_set elements will be considered, and you need to provide a
delta_t_minutes specification of the number of images and their cadence.

The final field, "microlensing_event_parameters", if provided, is used to configure
PSPL microlensing events to be injected into random stars. Note that these are
injected into the synthpop catalog, which may be smaller or larger than the output image size.

One output image is made for each epoch.

