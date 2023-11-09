"""
Input Parameters describing sims saved in ../photosphere/
"""

# imports and internal settings
import numpy as np
from astropy import units
Msun_str = units.Msun.to_string('latex_inline')


unitsOut = {
    'dist': units.au,
    'time': units.year,
    'temp': units.K,
    'density': units.g / units.cm**3,
}


# script input parameters
JOB_PROFILES_LIST = (
    {
        'job_name': '../photosphere/luis_2md/light',
        'file_indexes': np.arange(0, 17600+1, 100),
        'plot_title_suffix': f" for 1.7{Msun_str} primary star with Nucleation Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'X': 0.691,
        'name': f"Dusty 1.7{Msun_str}",
        'nickname': '4md',
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'job_name': '../photosphere/luis_4md/light',
        'file_indexes': np.arange(0, 17600+1, 100),
        'plot_title_suffix': f" for 3.7{Msun_str} primary star with Nucleation Dust",
        'ieos': 10,
        'params': {
            'X' : 0.686,
            'Z' : 0.024,
            'mu': 2.381,
        },
        'X': 0.686,
        'name': f"Dusty 3.7{Msun_str}",
        'nickname': '4md',
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'job_name': '../photosphere/miguel_2m/binary',
        'file_indexes': np.arange(0, 5000+1, 100),
        'plot_title_suffix': f" for 1.7{Msun_str} primary star without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'X': 0.691,
        'name': f"Non-dusty 1.7{Msun_str}",
        'nickname': '2m',
        #'fmt': '--',
        'color': 'orange',
    },
    {
        'job_name': '../photosphere/miguel_4m/binary',
        'file_indexes': np.arange(0, 5000+1, 100),
        'plot_title_suffix': f" for 3.7{Msun_str} primary star without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.686,
            'Z' : 0.024,
            'mu': 2.381,
        },
        'X' : 0.686,
        'name': f"Non-dusty 3.7{Msun_str}",
        'nickname': '4m',
        #'fmt': '--',
        'color': 'orange',
    },
    {
        'job_name': '../photosphere/miguel_2m_2022/binary',
        'file_indexes': np.arange(0, 6000+1, 100),
        'plot_title_suffix': f" for 1.7{Msun_str} primary star without Dust",
        'ieos': 10,
        'params': {
            'X' : 0.691,
            'Z' : 0.021,
            'mu': 2.381,
        },
        'X': 0.691,
        'name': f"Non-dusty 1.7{Msun_str}",
        'nickname': '2m_2022',
        #'fmt': '--',
        'color': 'orange',
    },
)

JOB_PROFILES_GROUPS = {
    '2m': (JOB_PROFILES_LIST[0], JOB_PROFILES_LIST[2],),
    '4m': (JOB_PROFILES_LIST[1], JOB_PROFILES_LIST[3],),
}
