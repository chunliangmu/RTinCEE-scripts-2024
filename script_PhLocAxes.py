#!/usr/bin/env python
# coding: utf-8

# # Main

# In[2]:


"""Scripts for analyzing of phantom outputs.

This script writes json files for each dump (and one json file synthsizing all outputs)
    to plot photosphere size vs time or orbital separation.
It does so by plotting photosphere intersection with traced rays originating from the primary star
    and shooting along the axes of the coordination frame.

"""


# ## Imports & Settings

# In[3]:


#%matplotlib inline
import math
import numpy as np
from numpy import pi
#import pandas
from astropy import units
from astropy import constants as const
import matplotlib.pyplot as plt
import matplotlib as mpl
#from moviepy.editor import ImageSequenceClip


# In[4]:


# import modules listed in ./lib/

from lib import clmuphantomlib as mupl
from lib.clmuphantomlib.readwrite import json_load, json_dump
from lib.clmuphantomlib.settings import DEFAULT_SETTINGS as settings
from lib.clmuphantomlib.log import error, warn, note, debug_info


#     ## import modules in arbitrary directory
#     
#     #import sys
#     
#     ## path to my python module lib directory
#     ## *** CHECK THIS! *** #
#     #SRC_LIB_PATH = sys.path[0] + '/lib'
#     
#     #if SRC_LIB_PATH not in sys.path:
#     #    sys.path.append(SRC_LIB_PATH)
#     ##print(*sys.path, sep='\n')    # debug
#     #print(
#     #    "\n*   Please Make sure my module files are located in this directory (or change the SRC_LIB_PATH variable):",
#     #    f"\n{SRC_LIB_PATH = }\n"
#     #)

# In[5]:


# parallels & optimizations


#import os
## Fixing stupid numba killing kernel
## See here https://github.com/numba/numba/issues/3016
#os.environ['NUMBA_DISABLE_INTEL_SVML']  = '1'
#from numba import njit, prange


from multiprocessing import cpu_count, Pool #Process, Queue

NPROCESSES = cpu_count()
if NPROCESSES is None:
    NPROCESSES = 1
NPROCESSES = max(NPROCESSES, 1)


# In[6]:


# settings
#
#   imported from script_input.py file

from script_PhLocAxes__input import iverbose, PHOTOSPHERE_TAU, JOB_PROFILES


# set metadata
with open("_metadata__input.json", 'r') as f:
    metadata = mupl.json_load(f)
metadata['Title'] = "Getting photosphere size on x, y, z axes"
metadata['Description'] = f"""Tracing 6 rays on +x, -x, +y, -y, +z, -z directon and get photosphere size, h, rho, u, T from them."""


plt.rcParams.update({'font.size': 20})


# print debug info
if iverbose >= 2:
    print(f"   Note: Will use {NPROCESSES} processes for parallelization")
    


# # Analysis

# ## Photosphere size vs time

# In[6]:


def write_ph_loc_axes(
    job_profile : dict,
    file_indexes : np.ndarray,
    rays_dir_def : dict,    # dict of list
    eos : mupl.eos_base.EoS_Base,
    photosphere_tau = PHOTOSPHERE_TAU,
    iverbose : int = 2,
):

    """Writing the photosphere locations of each dump to json files.

    Notes:
    Using mpdf.params['hfact']
    """
    
    
    mpdf = mupl.MyPhantomDataFrames()

    
    job_name = job_profile['job_name']
    X = job_profile['X']
    ieos = job_profile['ieos']

    
    # init rays directions
    rays_dir = {}
    for key in rays_dir_def.keys():
        rays_dir[key] = np.array(rays_dir_def[key])


    # main
    for file_index in file_indexes:
        
        # init answer dict / array
        photosphere_pars = { # [legend][par_name][time]
            'time_yr': None,
            'orbsep_Rsun': None,
            'data': {},
            'rays_dir': rays_dir_def,
            'rays': {},
        }  
        for key in rays_dir.keys():
            photosphere_pars['data'][key] = {}

        # read data
        mpdf.read(job_name, file_index, reset_xyz_by='CoM', iverbose=iverbose)
        if 'Tdust' in mpdf.data['gas'].columns:
            mpdf.data['gas']['T'] = mpdf.data['gas']['Tdust']
        elif 'temperature' in mpdf.data['gas'].columns:
            mpdf.data['gas']['T'] = mpdf.data['gas']['temperature']
        mpdf.calc_sdf_params(
            calc_params=['kappa',], #'R1',
            calc_params_params={'ieos': ieos, 'X':X, 'overwrite':False, 'kappa_translate_from_cgs_units':True},
            iverbose=iverbose,
        )
        hfact = mpdf.params['hfact']
        mpart = mpdf.params['mass']
        
        photosphere_pars['time_yr'] = mpdf.get_time().to_value(units.year)
        photosphere_pars['orbsep_Rsun'] = mpdf.get_orb_sep().to_value(units.Rsun)

        # construct rays_dict
        star_loc = np.array([mpdf.data['sink'][axis][0] for axis in 'xyz'])
        rays_dict = {}    # legend: ray
        for key in rays_dir.keys():
            # init
            ray = np.array([
                star_loc,
                star_loc + rays_dir[key],
            ])
            rays_dict[key] = ray
            photosphere_pars['rays'][key] = ray.tolist()
            ray_unit_vec = ray[1, :] - ray[0, :]
            ray_unit_vec = ray_unit_vec / np.sum(ray_unit_vec**2)**0.5


            # optimization- first select only the particles affecting the ray
            #  because interpolation of m points with N particles scales with O(N*m),
            #  reducing N can speed up calc significantly
            sdf = mpdf.data['gas']
            kernel_radius = sdf.kernel.get_radius()
            hs = np.array(sdf['h'])
            pts = np.array(sdf[['x', 'y', 'z']])    # (npart, 3)-shaped array
            pts_on_ray = mupl.get_closest_pt_on_line(pts, ray)
            sdf_selected_indices = (np.sum((pts - pts_on_ray)**2, axis=-1) <= (kernel_radius * hs)**2)
            if iverbose:
                debug_info(
                    'write_ph_loc_axes()', iverbose,
                    f"{np.count_nonzero(sdf_selected_indices)} particles are close enough to the ray to have effects."
                )
            sdf = sdf.iloc[sdf_selected_indices]
            pts = np.array(sdf[['x', 'y', 'z']])    # (npart, 3)-shaped array


            # get optical depth
            if iverbose:
                debug_info(
                    'write_ph_loc_axes()', iverbose,
                    f"{ray = }"
                )
            pts_on_ray, dtaus, pts_order = mupl.get_optical_depth_by_ray_tracing_3D(sdf=sdf, ray=ray)
            photosphere, waypts_list = mupl.get_photosphere_on_ray(
                pts_on_ray, dtaus, pts_order, sdf, ray,
                calc_params = ['loc', 'R1', 'rho', 'u', 'h', 'T'],
                hfact = hfact, mpart=mpart, eos=eos, sdf_units=mpdf.units,
                ray_unit_vec=ray_unit_vec,
            )
            photosphere_pars['data'][key] = photosphere
            photosphere_pars['data'][key]['size'] = photosphere['R1']
                
            if iverbose:
                debug_info(    # debug
                    'write_ph_loc_axes()', iverbose,
                    f"{photosphere_loc = }\n{photosphere_dist_to_ray0 = }\n",
                    f"{photosphere_taus = }\n",
                    f"{pts_on_ray_ordered[photosphere_loc_index:photosphere_loc_index+2] = }",
                )

        with open(f"{mpdf.get_filename()}__photospherePars__xyz.json", 'w') as f:
            json_dump(photosphere_pars, f, metadata=metadata)
            if iverbose: print(f"\n\nWritten to {f.name}\n")

    return None


# In[9]:


# main process



# init rays directions
rays_dir_def = {
    # legend: ray direction name
    '+x'  : [ 1., 0., 0.],
    '+y'  : [ 0., 1., 0.],
    '+z'  : [ 0., 0., 1.],
    '-x'  : [-1., 0., 0.],
    '-y'  : [ 0.,-1., 0.],
    '-z'  : [ 0., 0.,-1.],
}


# run main

if __name__ == '__main__':
    
    
    # get ph loc for each dump file
    for job_profile in JOB_PROFILES:
    
        file_indexes = job_profile['file_indexes']
        eos = mupl.get_eos(job_profile['ieos'], job_profile['params'], settings)
    
        
        if NPROCESSES <= 1:
            
            # single process
    
            write_ph_loc_axes(
                job_profile = job_profile, file_indexes = file_indexes, rays_dir_def = rays_dir_def,
                eos = eos, photosphere_tau = PHOTOSPHERE_TAU, iverbose = iverbose,
            )
            
        else:
            
            # multi-process

            args = [
                (job_profile,
                file_index,
                rays_dir_def,
                eos,
                PHOTOSPHERE_TAU,
                0,
                ) for file_index in file_indexes
            ]

            with Pool(processes=NPROCESSES) as pool:
                pool.starmap(write_ph_loc_axes, args)
    
    


# In[12]:


file_indexes = JOB_PROFILES[0]['file_indexes']
file_indexes[:, np.newaxis].tolist()


# In[32]:


if __name__ == '__main__':
    
    # syntheize the files into one big file
    
    for job_profile in JOB_PROFILES:
    
        job_name     = job_profile['job_name']
        file_indexes = job_profile['file_indexes']
    
    
        # init
        photosphere_pars_all = { # [legend][par_name][time]
            'time_yr': [],
            'orbsep_Rsun': [],
            'data': {},
            'rays_dir': rays_dir_def,
            'rays': {},
        }  
        for key in rays_dir_def.keys():
            photosphere_pars_all['data'][key] = {
                'size': [],
                'rho' : [],
                'u'   : [],
                'h'   : [],
                'T'   : [],
            }
            photosphere_pars_all['rays'][key] = []
    
        
        # fetch
        for file_index in file_indexes:
            with open(f"{job_name}_{file_index:05}__photospherePars__xyz.json", 'r') as f:
                
                if iverbose: print(f"\n\nLoading {f.name}... ", end='')
                
                photosphere_pars = json_load(f)
                for it in ['time_yr', 'orbsep_Rsun']:
                    photosphere_pars_all[it].append(photosphere_pars[it])
                for key in rays_dir_def.keys():
                    for it in photosphere_pars_all['data'][key].keys():
                        obj = photosphere_pars['data'][key][it]
                        if isinstance(obj, np.ndarray):
                            obj = obj.tolist()
                        photosphere_pars_all['data'][key][it].append(obj)
                    photosphere_pars_all['rays'][key].append(photosphere_pars['rays'][key]) 
    
                if iverbose: print(f"Done.\n")
    
        
        # write
        with open(f"{job_name}__photospherePars__xyz.json", 'w') as f:
            json_dump(photosphere_pars_all, f, metadata=metadata)
            if iverbose: print(f"\n\nWritten to {f.name}.\n")


    print("\n\n\n*** All Done. ***\n\n\n")


# In[ ]:




