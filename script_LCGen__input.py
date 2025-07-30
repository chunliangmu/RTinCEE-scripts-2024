"""
_photosphere_jobProfiles__inputInput Parameters for my scripts.
"""

# imports and internal settings
import numpy as np
from numpy import pi
from astropy import units
from astropy import constants as const
from _photosphere_jobProfiles__input import interm_dir, output_dir, unitsOut, SPEC_DIST, PHOTOSPHERE_TAU, AVG_KC_PP, JOB_PROFILES_DICT #, fps

unitsOut['flux_wav'] = (units.erg / units.s / units.cm**2) / units.angstrom


verbose = 6

job_nicknames = ['2mdd', '2mdc', '2mds', '2mdh',]#, '2md', '4m', '2m_2022', '2md']
xyzs_list  = ['xyz', 'xzy', 'yzx']
no_xy=(256, 256)
no_xy_txt = 'x'.join([f'{i}' for i in no_xy])
output_dir = f'../fig/20240222_LCGen/{no_xy_txt}/'
interm_dir += 'olim_'
verbose_loop = 0

nsample_pp  : int   = 1000   # no of sample points per particle for integration
z_olim_kc   : float = 1.152  # col kernel limit for when srcfunc began to count for olim

# at t=0... (only used when use_Tscales=True)
# numbers from Gonzalez-2024-1
use_Tscales : None|str = None #None 'scale', 'cut', 'delete'
if use_Tscales: interm_dir += f'T{use_Tscales}_' #'Tscaled_'
Ls_mesa = {    # from Gonzalez-2024-1
    '2md': 5180 * units.Lsun,
    '4md': 1.19e4 * units.Lsun,
}

Rs_ph_mesa = {    # from Gonzalez-2024-1
    '2mdnr': 260 * units.Rsun,
    '2md': 260 * units.Rsun,
    '4md': 330 * units.Rsun,
}

if use_Tscales == 'scale':
    # Set to make initial luminosity close to mesa value
    #   obtained by iterating Rs_ph_init at t=0 until L ~ L_mesa
    Rs_ph_init = {
        '2md': 359 * units.Rsun,
        '4md': 464 * units.Rsun,
    }
elif use_Tscales == 'cut':
    Rs_ph_init = {
        '2md': 292 * units.Rsun,
        '4md': 396 * units.Rsun,
    }
else:
    Rs_ph_init = Rs_ph_mesa

Ts_ph_init = {
    n: ((Ls_mesa[n] / (4*pi*Rs_ph_init[n]**2 * const.sigma_sb))**(1/4)).to(units.K)
    for n in Ls_mesa.keys()
}
Ts_ph_init['2mdnr'] = 3227*units.K

# - SED settings -
# freq: minimum range 1e9~1e20 Hz (covering microwave to x-ray)
# wavelen 1e-11m ~ 0.1m
#wavlens = (np.logspace(-10, 0, 10000) * units.m).cgs   # default (broad)
wavlens = (np.logspace(-2, 5., 10000) * units.um).cgs   # default
#wavlens = (np.logspace(-1, 3.5, 50) * units.um).cgs  # mcfost


if __name__ == '__main__':
    # gen Tscales
    from _sharedFuncs import gen_Tscales
    if verbose: print('note', None, verbose, f"\n{Ts_ph_init = }\n")
    for job_nickname in job_nicknames:
        job_profile = JOB_PROFILES_DICT[job_nickname]
        scales = gen_Tscales(
            job_name = job_profile['job_name'],
            T_ph = Ts_ph_init[job_nickname],
            R_ph = Rs_ph_mesa[job_nickname],
            method = use_Tscales,
            do_save= True,
            params = job_profile['params'],
            verbose=verbose,
        )
        
