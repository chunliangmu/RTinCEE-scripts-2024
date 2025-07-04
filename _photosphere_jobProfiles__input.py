#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Input Parameters describing sims saved in ../photosphere/
"""

# imports and internal settings
import numpy as np
from astropy import units
Msun_str = units.Msun.to_string('latex_inline')

# output directory- remember to add the final sep '/' at the end
interm_dir = '../interm/'
output_dir = '../fig/'
out_dir    = output_dir  # deprecated keyword


fps = 6

unitsOut = {
    'dist': units.au,
    'time': units.year,
    'temp': units.K,
    'density': units.g / units.cm**3,
    'opacity': units.cm**2 / units.g,
    'dimless': units.dimensionless_unscaled,
    'speed': units.km/units.s,
    'flux' : (units.erg / units.s / units.cm**2),
    #'flux_wav' : (units.erg / units.s / units.cm**2) / units.angstrom,
}

SPEC_DIST = 10 * units.parsec
PHOTOSPHERE_TAU = 2./3. #np.log(2)


# script input parameters
JOB_PROFILES_LIST = (
    {
        'raw_dir' : '../raw/clmu_2mdh/',
        'file_prefix': 'lumo',
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star (halfway companion, forced MESA profile)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Halfway 1.7{Msun_str}",
        'nickname': '2mdh',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mds/',
        'file_prefix': 'lumo',
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star (surface companion, forced MESA profile)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Surface 1.7{Msun_str}",
        'nickname': '2mds',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdd/',
        'file_prefix': 'lumo',
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star (donor only, forced MESA profile)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Donor 1.7{Msun_str}",
        'nickname': '2mdd',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdc/',
        'file_prefix': 'lumo',
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star with forced MESA profile",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Cut 1.7{Msun_str}",
        'nickname': '2mdc',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_sol/',
        'file_prefix': 'sol',
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for Sun without full relaxation",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 0.5988,
        },
        'name': f"NoRelax Sun",
        'nickname': 'sol',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnr/',
        'file_prefix': '2mdnr', # nr means no-relax
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star without full relaxation",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str}",
        'nickname': '2mdnr',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnrs/',
        'file_prefix': '2mdnrs', # s for surface
        'file_indexes': np.arange(10),
        'plot_title_suffix': f" for 1.7{Msun_str} star without full relaxation, with companion star at the donor surface",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str} Surface",
        'nickname': '2mdnrs',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnrh/',
        'file_prefix': '2mdnrh', # h for halfway
        'file_indexes': np.arange(10),
        'plot_title_suffix': f" for 1.7{Msun_str} star without full relaxation, with companion star halfway between donor surface and RLOF",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str} Halfway",
        'nickname': '2mdnrh',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnrt0e1/',
        'file_prefix': '2mdnrt0e1', # t=0, exp 1
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star without full relaxation res x8 (npart=11M)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str} x8 npart",
        'nickname': '2mdnrt0e1',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnrt0e2/',
        'file_prefix': '2mdnrt0e2', # t=0, exp 2
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star without full relaxation res x64 (npart=88M)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str} x64 npart",
        'nickname': '2mdnrt0e2',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/clmu_2mdnrt0e2n/',
        'file_prefix': '2mdnrt0e2n', # t=0, exp 2, no relax-o-matic
        'file_indexes': np.array([0]),
        'plot_title_suffix': f" for 1.7{Msun_str} star without relaxation res x64 (npart=88M)",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"NoRelax 1.7{Msun_str} x64 npart",
        'nickname': '2mdnrt0e2n',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/luis_2md/',
        'file_prefix': 'light',
        #'job_name': '../photosphere/luis_2md/light', # deprecated keyword- will still be added automatically later
        #'file_indexes': np.arange(0, 17600+1, 100),
        'file_indexes': np.concatenate((np.arange(0, 5000-1, 20), np.arange(5000, 17600+1, 50))),
        'plot_title_suffix': f" for 1.7{Msun_str} primary with Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Dusty 1.7{Msun_str}",
        'nickname': '2md',
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/luis_4md/',
        'file_prefix': 'light',
        #'job_name': '../photosphere/luis_4md/light',
        'file_indexes': np.arange(0, 17600+1, 100),
        'plot_title_suffix': f" for 3.7{Msun_str} primary with Dust",
        'ieos': 10,
        'params': {
            'X' : 0.686,
            'Z' : 0.024,
            'mu': 2.381,
        },
        'name': f"Dusty 3.7{Msun_str}",
        'nickname': '4md',
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'raw_dir' : '../raw/miguel_2m/',
        'file_prefix': 'binary',
        #'job_name': '../photosphere/miguel_2m/binary',
        'file_indexes': np.arange(0, 5000+1, 100),
        'plot_title_suffix': f" for 1.7{Msun_str} primary without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Non-dusty 1.7{Msun_str}",
        'nickname': '2m',
        #'fmt': '--',
        'color': 'orange',
    },
    {
        'raw_dir' : '../raw/miguel_4m/',
        'file_prefix': 'binary',
        #'job_name': '../photosphere/miguel_4m/binary',
        'file_indexes': np.arange(0, 5000+1, 100),
        'plot_title_suffix': f" for 3.7{Msun_str} primary without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.686,
            'Z' : 0.024,
            'mu': 2.381,
        },
        'name': f"Non-dusty 3.7{Msun_str}",
        'nickname': '4m',
        #'fmt': '--',
        'color': 'orange',
    },
    {
        'raw_dir' : '../raw/miguel_2m_2022/',
        'file_prefix': 'binary',
        #'job_name': '../photosphere/miguel_2m_2022/binary',
        'file_indexes': np.arange(0, 6000+1, 100),
        'plot_title_suffix': f" for 1.7{Msun_str} primary without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'name': f"Non-dusty 1.7{Msun_str}",
        'nickname': '2m_2022',
        #'fmt': '--',
        'color': 'orange',
    },
)

for job_profile in JOB_PROFILES_LIST:
    # deprecated keyword generation
    job_profile['job_name'] = job_profile['raw_dir'] + job_profile['file_prefix']
    job_profile['X'] = job_profile['params']['X']

JOB_PROFILES_GROUPS = {
    '2m': (JOB_PROFILES_LIST[0], JOB_PROFILES_LIST[2],),
    '4m': (JOB_PROFILES_LIST[1], JOB_PROFILES_LIST[3],),
}

JOB_PROFILES_DICT = { job_profile['nickname']: job_profile for job_profile in JOB_PROFILES_LIST }
