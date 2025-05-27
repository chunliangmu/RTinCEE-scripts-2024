"""
Input Parameters for my scripts.
"""

# imports and internal settings
from _photosphere_jobProfiles__input import (
    output_dir, interm_dir,
    fps, unitsOut,
    PHOTOSPHERE_TAU,
    JOB_PROFILES_LIST, JOB_PROFILES_GROUPS, JOB_PROFILES_DICT,
)
from clmuphantomlib.units_util import DEFAULT_UNITS


# script input parameters

verbose = 3

interm_dir += 'test_'
output_dir = f'../fig/20241212_photosphere-axes_v-vs-R1/test_'


# units used in the dump file
# *** WARNING: script does not check this yet!
#      Change it if you are not using Rsun & Msun for phantom dumps!!!
unitsIn = DEFAULT_UNITS

JOB_PROFILES = JOB_PROFILES_LIST[:]
