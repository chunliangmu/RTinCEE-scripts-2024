"""
Input Parameters for my scripts.
"""

# imports and internal settings
import numpy as np
from astropy import units
from _photosphere_jobProfiles__input import interm_dir, output_dir, unitsOut, SPEC_DIST, PHOTOSPHERE_TAU, JOB_PROFILES_DICT #, fps

unitsOut['flux_wav'] = (units.erg / units.s / units.cm**2) / units.angstrom


verbose = 6

job_nicknames = ['2md', '4md']#, '4m', '2m_2022']
xyzs_list  = ['xyz', 'xzy']
no_xy=(64, 64)
no_xy_txt = 'x'.join([f'{i}' for i in no_xy])
output_dir = f'../fig/20240222_LCGen/{no_xy_txt}/test_'
verbose_loop = 0

# - SED settings -
# freq: minimum range 1e9~1e20 Hz (covering microwave to x-ray)
# wavelen 1e-11m ~ 0.1m
#wavlens = (np.logspace(-10, 0, 10000) * units.m).cgs   # default (broad)
wavlens = (np.logspace(-2, 5., 10000) * units.um).cgs   # default
#wavlens = (np.logspace(-1, 3.5, 50) * units.um).cgs  # mcfost