#!/usr/bin/env python
# coding: utf-8

# # Main

# In[7]:


"""Scripts for analyzing of phantom outputs.

This script writes json files for each dump (and one json file synthsizing all outputs)
    to plot photosphere size vs time or orbital separation.
It does so by plotting photosphere intersection with traced rays originating from the primary star
    and shooting along the axes of the coordination frame.

"""


# ## Imports & Settings

# In[8]:


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
import json


# In[9]:


# import modules listed in ./lib/

from lib import clmuphantomlib as mupl


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

# In[4]:


# parallels & optimizations


#import os
## Fixing stupid numba killing kernel
## See here https://github.com/numba/numba/issues/3016
#os.environ['NUMBA_DISABLE_INTEL_SVML']  = '1'
#from numba import njit, prange


from multiprocessing import cpu_count, Process, Queue

NPROCESSES = cpu_count()
if NPROCESSES is None:
    NPROCESSES = 1
NPROCESSES = max(NPROCESSES, 1)


# In[5]:


# settings
#
#   imported from script_input.py file

from script_PhLocAxes__input import iverbose, PHOTOSPHERE_TAU, JOB_PROFILES


plt.rcParams.update({'font.size': 20})


# print debug info
if iverbose >= 2:
    print(f"   Note: Will use {NPROCESSES} processes for parallelization")
    


#     # Test
# 
#     photosphere_tau = PHOTOSPHERE_TAU
#     iverbose = 3
# 
#     job_profile = JOB_PROFILES[0]
#     job_name = job_profile['job_name']
#     ieos = job_profile['ieos']
#     file_indexes = job_profile['file_indexes']
#     plot_title_suffix = job_profile['plot_title_suffix']
#     X = job_profile['X']
# 
# 
#     file_index=file_indexes[5]
#     mpdf = mupl.MyPhantomDataFrames()
#     mpdf.read(
#         job_name, file_index,
#         calc_params=['T', 'kappa', 'R1'],
#         reset_xyz_by="R1",
#         calc_params_params={'ieos': ieos, 'X':X, 'overwrite':False, 'kappa_translate_from_cgs_units':True},
#         iverbose=iverbose,
#     )
#     mpdf.plot_render(plot_title_suffix=plot_title_suffix,
#         xlim=(-60000, 60000), ylim=(-60000, 60000),
#         norm=mpl.colors.LogNorm(vmin=1e-25, vmax=1e-5, clip=True),
#     )
#     if iverbose:
#         print()
#         print(mpdf.get_time())
#         print(mpdf.data['gas'].keys())
#         print(mpdf.data['sink'])

# 
#     
#     
#     # units info (for plot axes title)
#     
#     unitsIn = {    # not used for reading phantom data dumps since it records units in the dump already
#         'dist': units.solRad,
#         'temp': units.K,
#         'flux': units.erg/units.cm**2/units.s,    # equivalent of units.g/units.s**3
#         'opacity': units.cm**2/units.g,
#     }
#     
#     unitsOut = {
#         'dist': units.solRad,
#         'temp': units.K,
#         'flux': units.erg/units.cm**2/units.s,    # equivalent of units.g/units.s**3
#         'opacity': units.cm**2/units.g,
#     }
#     
#     unitsOutTxt = {}
#     
#     for key in unitsOut.keys():
#         unitsOutTxt[key] = unitsOut[key].to_string('latex_inline')
# 

# # Analysis

# ## Photosphere size vs time

# In[6]:


def write_ph_loc_axes(
    job_profile : dict,
    file_indexes : np.ndarray,
    rays_dir_def : dict,    # dict of list
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
            photosphere_pars['data'][key] = {
                'size': None,
                'rho' : None,
                'u'   : None,
                'h'   : None,
                'T'   : None,
            }

        # read data
        mpdf.read(job_name, file_index, reset_xyz_by='CoM', iverbose=iverbose)
        if 'Tdust' in mpdf.data['gas'].columns:
            mpdf.data['gas']['T'] = mpdf.data['gas']['Tdust']
        elif 'temperature' in mpdf.data['gas'].columns:
            mpdf.data['gas']['T'] = mpdf.data['gas']['temperature']
        mpdf.calc_sdf_params(
            calc_params=['T', 'kappa',], #'R1',
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
            if iverbose >= 3:
                print(f"    Info: {np.count_nonzero(sdf_selected_indices)} particles are close enough to the ray to have effects.")
            sdf = sdf.iloc[sdf_selected_indices]
            pts = np.array(sdf[['x', 'y', 'z']])    # (npart, 3)-shaped array


            # get optical depth
            if iverbose >= 3: print(f"    Info: {ray = }")
            pts_on_ray, dtaus, pts_order = mupl.get_optical_depth_by_ray_tracing_3D(sdf=sdf, ray=ray)
            pts_order_nonzero = np.where(dtaus[pts_order])[0]
            if iverbose >= 3:
                print(f"    Info: {pts_order_nonzero.size = }")
            pts_on_ray_ordered = pts_on_ray[pts_order]
            dist_to_ray0_ordered = np.sum((pts_on_ray_ordered - ray[0]) * ray_unit_vec, axis=-1)
            taus_ordered = np.cumsum(dtaus[pts_order])


            # get photosphere
            photosphere_loc_index = np.searchsorted(taus_ordered, photosphere_tau) - 1
            photosphere_found = photosphere_loc_index <= len(sdf) - 2
            if photosphere_found:
                if photosphere_loc_index == -1:
                    # if first particle blocks everything: estimate it to be where that particle is
                    photosphere_loc_index = 0
                    photosphere_loc = pts_on_ray_ordered[0]
                    photosphere_taus = [0., taus_ordered[photosphere_loc_index]]    # for debug
                else:
                    # intepolate to find loc
                    photosphere_taus = taus_ordered[photosphere_loc_index : photosphere_loc_index+2]
                    photosphere_dtau = photosphere_taus[1] - photosphere_taus[0]
                    photosphere_dtau0_frac = (photosphere_taus[1] - photosphere_tau) / photosphere_dtau
                    photosphere_loc = \
                        pts_on_ray_ordered[photosphere_loc_index] * photosphere_dtau0_frac + \
                        pts_on_ray_ordered[photosphere_loc_index+1] * (1 - photosphere_dtau0_frac)
                photosphere_dist_to_ray0 = np.sum((photosphere_loc - ray[0]) * ray_unit_vec, axis=-1)


                photosphere_rho = mupl.get_sph_interp(sdf, 'rho' , photosphere_loc)
                photosphere_u   = mupl.get_sph_interp(sdf, 'u'   , photosphere_loc)
                
                photosphere_pars['data'][key]['size'] = photosphere_dist_to_ray0
                photosphere_pars['data'][key]['rho' ] = photosphere_rho
                photosphere_pars['data'][key]['u'   ] = photosphere_u
                photosphere_pars['data'][key]['h'   ] = hfact * (mpart / photosphere_rho)**(1./3.)
                photosphere_pars['data'][key]['T'   ] = mupl.get_temp_from_u(
                    rho=photosphere_rho, u=photosphere_u, mu=None, ieos=ieos,
                    rho_unit=mpdf.units['density'], u_unit=mpdf.units['specificEnergy'],
                ).item()
                
                if iverbose >= 3:
                    print(    # debug
                        f"{photosphere_loc = }\n{photosphere_dist_to_ray0 = }\n",
                        f"{photosphere_taus = }\n",
                        f"{pts_on_ray_ordered[photosphere_loc_index:photosphere_loc_index+2] = }",
                    )
            else:
                if iverbose: print(f"*    Warning: Photosphere not found ({photosphere_loc_index=}; {max(taus_ordered)=})")
            
        with open(f"{mpdf.get_filename()}__photospherePars__xyz.json", 'w') as f:
            json.dump(photosphere_pars, f)
            if iverbose: print(f"\n\nWritten to {f.name}\n")

    return None


# In[7]:


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
    
        
        if NPROCESSES <= 1:
            
            # single process
    
            write_ph_loc_axes(
                job_profile = job_profile, file_indexes = file_indexes, rays_dir_def = rays_dir_def,
                photosphere_tau = PHOTOSPHERE_TAU, iverbose = 0,
            )
            
        else:
            
            # multi-process
            
            file_indexes_list = np.array_split(file_indexes, NPROCESSES)
            processes_list = []
            for i, file_indexes in enumerate(file_indexes_list):
                processes_list.append(
                    Process(
                        target=write_ph_loc_axes,
                        kwargs={'job_profile': job_profile, 'file_indexes': file_indexes, 'rays_dir_def': rays_dir_def,
                                'photosphere_tau': PHOTOSPHERE_TAU, 'iverbose': 0,
                        },
                    )
                )
        
            for process in processes_list:
                process.start()
        
            for process in processes_list:
                process.join()
    
    


# In[8]:


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
                
                photosphere_pars = json.load(f)
                for it in ['time_yr', 'orbsep_Rsun']:
                    photosphere_pars_all[it].append(photosphere_pars[it])
                for key in rays_dir_def.keys():
                    for it in photosphere_pars_all['data'][key].keys():
                        photosphere_pars_all['data'][key][it].append(photosphere_pars['data'][key][it])
                    photosphere_pars_all['rays'][key].append(photosphere_pars['rays'][key]) 
    
                if iverbose: print(f"Done.\n")
    
        
        # write
        with open(f"{job_name}__photospherePars__xyz.json", 'w') as f:
            json.dump(photosphere_pars_all, f)
            if iverbose: print(f"\n\nWritten to {f.name}.\n")

