"""
Input Parameters for my scripts.
"""

# imports and internal settings
from _photosphere_jobProfiles__input import unitsOut, JOB_PROFILES_LIST, JOB_PROFILES_GROUPS


# script input parameters

verbose = 3

PHOTOSPHERE_TAU = 1.0

ray_no = 600

plane_axes_list = ['xyz', 'xzy']

box_lim = 199
box_lim_dict = { # using JOB_PROFILES_LIST['nickname'] as key
    '2md': 199,
    '4md': 199,
    '2m' : 39,
    '4m' : 39,
    '2m_2022': 39,
}

fps = 10

use_saved_jsons = False

JOB_PROFILES = JOB_PROFILES_LIST[:]
