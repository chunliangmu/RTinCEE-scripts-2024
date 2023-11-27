#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""Scripts for analyzing of phantom outputs.

This script analyze the dump files about the opacity based on density and temperature.

"""


# # Main

# ## Imports & Settings

# In[2]:


import numpy as np
from astropy import units
import matplotlib.pyplot as plt
import matplotlib as mpl
from moviepy.editor import ImageSequenceClip
from os import path


# In[3]:


# import my modules listed in ./main/

from main import clmuphantomlib as mupl
#from main.clmuphantomlib.readwrite import json_load
from main.clmuphantomlib.settings   import DEFAULT_SETTINGS as settings
from main.clmuphantomlib.units_util import get_val_in_unit #set_as_quantity, get_units_field_name, get_units_cgs
from main.clmuphantomlib.eos_mesa   import EoS_MESA_opacity


# In[4]:


# settings
#
#   imported from script_input.py file


from script_kappaProfile__input import verbose, fps, unitsOut, JOB_PROFILES_DICT

unitsOutTxt = {  key  : unitsOut[key].to_string('latex_inline') for key in unitsOut.keys() }


plt.rcParams.update({'font.size': 20})



# In[ ]:


# plotting kappa vs Temp and kappa vs rho
if __name__ == '__main__':

    ylim = (1e-6, 1e4)
    xlims= {
        'T'  : (  5.,  2e6),
        'rho': (2e-20, 1e-3),
    }
    
    mpdf = mupl.MyPhantomDataFrames()
    for key in ['2md', '4md']:
        job_profile = JOB_PROFILES_DICT[key]
        job_name    = job_profile['job_name']
        params      = job_profile['params']
        eos_opacity = EoS_MESA_opacity(params, settings)

        outfilenames = []
        for file_index in job_profile['file_indexes']:
            fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
            fig.subplots_adjust(wspace=0.0)
            mpdf.read(job_name, file_index, reset_xyz_by_CoM=True, verbose=verbose)
            #jobfilename = mupl.get_filename_phantom_dumps(job_name, file_index)
            jobfilename = mpdf.get_filename()
            # get temperature column label (one of the elem in the set below)
            temp_key = {'T', 'temperature', 'Tdust'}.intersection(mpdf.data['gas'].keys()).pop()
            mpdf.data['gas']['T'    ] = mpdf.data['gas'][temp_key]
            mpdf.data['gas']['kappa'] = get_val_in_unit(mpdf.data['gas']['kappa'], units.cm**2/units.g, mpdf.units['opacity'])
    
            y_orig = mpdf.get_val('kappa').to(unitsOut['opacity'])
            y_mesa = eos_opacity.get_kappa(mpdf.get_val('rho'), mpdf.get_val('T'), do_extrap=False).to(unitsOut['opacity'])
            y_extrap_indexes = np.where(~np.isfinite(y_mesa))[0]
            y_mesa_extrap = eos_opacity.get_kappa(
                mpdf.get_val('rho')[y_extrap_indexes], mpdf.get_val('T')[y_extrap_indexes], do_extrap=True).cgs
    
            
            # kappa vs temp
            ax = axes[0]
            x = mpdf.get_val('T').to(unitsOut['temp'], equivalencies=units.equivalencies.temperature())
            x_extrap = x[y_extrap_indexes]
            ax.loglog(x, y_orig, '.')
            ax.loglog(x, y_mesa, '.')
            ax.loglog(x_extrap, y_mesa_extrap, '.')
            ax.set_xlim(xlims['T'])
            ax.set_ylim(ylim)
            ax.set_xlabel(f"$T$ / {x.unit.to_string('latex_inline')}")
            ax.set_ylabel(f"$\\kappa$ / {y_orig.unit.to_string('latex_inline')}")
            ax.text(
                0.02, 0.98,
                f"Time = {mpdf.get_time(unitsOut['time']):.1f}",
                color = "black", ha = 'left', va = 'top',
                transform=ax.transAxes,
            )
    
            # kappa vs rho
            ax = axes[1]
            x = mpdf.get_val('rho').to(unitsOut['density'])
            x_extrap = x[y_extrap_indexes]
            ax.loglog(x, y_orig, '.', label='Nucleation')
            ax.loglog(x, y_mesa, '.', label='MESA')
            ax.loglog(x_extrap, y_mesa_extrap, '.', label='MESA extrap')
            ax.set_xlim(xlims['rho'])
            ax.set_ylim(ylim)
            ax.set_xlabel(f"$\\rho$ / {x.unit.to_string('latex_inline')}")
            ax.legend(loc='lower right')
    
    
            fig.suptitle(
                f"Opacity of all particles in the dump (different calculation method)\n"
                #f"resolution = {mpdf.params['nparttot']:.2e}\n"
                f"{job_profile['plot_title_suffix']}"
            )
    
            outfilename = f"{jobfilename}__kappaProfile.png"
            fig.savefig(outfilename)
            outfilenames.append(outfilename)
            plt.close(fig)


        # define job_folder_prefix
        for i in range(len(job_name)-1, -1, -1):
            if job_name[i] == path.sep:
                job_folder_prefix = job_name[:i]
                break
            else:
                job_folder_prefix = job_name
        with ImageSequenceClip(outfilenames, fps=fps) as vid:
            moviefilename = f"{job_folder_prefix}__kappaProfile__movie.mp4"
            vid.write_videofile(moviefilename)


# if __name__ == '__main__':
#     
#     plt.close(fig)
#     
#     fig, ax = plt.subplots(figsize=(8, 8))
#     x = mpdf.get_val('T').to(unitsOut['temp'], equivalencies=units.equivalencies.temperature())
#     y = mpdf.get_val('rho').to(unitsOut['density'])
#     ax.loglog(x, y, '.')
#     ax.text(
#         0.95, 0.5,
#         f"Time = {mpdf.get_time(unitsOut['time']):.1f}",
#         color = "black", ha = 'right', va = 'top',
#         transform=ax.transAxes,
#     )
#     ax.set_xlim(xlims['T'])
#     ax.set_ylim(xlims['rho'])
# 
# 
# 
# 

# In[ ]:




