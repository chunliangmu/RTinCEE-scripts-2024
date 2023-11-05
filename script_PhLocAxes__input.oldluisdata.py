"""
Input Parameters for my scripts.
"""

# imports and internal settings
import numpy as np
from astropy import units
Msun_str = units.Msun.to_string('latex_inline')


# script input parameters

iverbose = 3

PHOTOSPHERE_TAU = 1.0

JOB_PROFILES_LIST = (
    {
        'job_name': '../photosphere/luis_2md/binary',
        'file_indexes': (0, 396, 792, 1188, 1584, 1980, 2178, 2376, 2574, 2772, 3564, 4710, 5780, 7090, 8000,),
        'plot_title_suffix': f" for 1.7{Msun_str} primary star with Nucleation Dust",
        'ieos': 10,
        'X': 0.691,
        'name': f"Dusty 1.7{Msun_str}",
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'job_name': '../photosphere/luis_4md/binary',
        'file_indexes': (
                0, 396, 792, 1188, 1584, 1980, 2178, 2277, 2376, 2574, 2772, 3564, 4710, 5780, 7090, 8000,
                11880, 15840, 17820, 19800,
            ),
        'plot_title_suffix': f" for 3.7{Msun_str} primary star with Nucleation Dust",
        'ieos': 10,
        'X': 0.686,
        'name': f"Dusty 3.7{Msun_str}",
        #'fmt': '-',
        'color': 'blue',
    },
    {
        'job_name': '../photosphere/miguel_2m/binary',
        'file_indexes': (0, 396, 792, 1187, 1584, 1978, 2376, 2772, 3564, 4709,),
        'plot_title_suffix': f" for 1.7{Msun_str} primary star without Dust",
        'ieos': 10,
        'X': 0.691,
        'name': f"Non-dusty 1.7{Msun_str}",
        #'fmt': '--',
        'color': 'orange',
    },
    {
        'job_name': '../photosphere/miguel_4m/binary',
        'file_indexes': (0, 396, 792, 1187, 1584, 1978, 2376, 2574, 2772, 3168, 3564, 4709,),
        'plot_title_suffix': f" for 3.7{Msun_str} primary star without Dust",
        'ieos': 10,
        'X': 0.686,
        'name': f"Non-dusty 3.7{Msun_str}",
        #'fmt': '--',
        'color': 'orange',
    },
)

JOB_PROFILES = JOB_PROFILES_LIST[0:2]