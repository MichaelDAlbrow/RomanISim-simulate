# RomanISim-simulate

Simulate images for the Nancy Grace Roman Space Telescope Time-Domain Survey.

This code uses synthpop and RomanISIM/STPSF to generate time-series simulated
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

```json
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

```

This is more-or-less the same as a standard synthpop config file, except only the
first l_set and b_set elements will be considered, and you need to provide a
delta_t_minutes specification of the number of images and their cadence.

The final field, "microlensing_event_parameters", if provided, is used to configure
PSPL microlensing events to be injected into random stars. Note that these are
injected into the synthpop catalog, which may be smaller or larger than the output image size.

One output image is made for each epoch.

