#!/usr/bin/env python
# coding: utf-8

"""Functions def for script_PhLocAxes

For some reason python native multiprocessing code throw erratic errors if functions not put in a separate file
    (See <https://stackoverflow.com/a/42383397>)
This aim to fix that.
"""

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

# import modules listed in ./lib/

#  import (my libs)
import clmuphantomlib as mupl
from clmuphantomlib import get_col_kernel_funcs
from clmuphantomlib.io import json_load, json_dump, hdf5_load, hdf5_dump
from clmuphantomlib.settings import DEFAULT_SETTINGS as settings
from clmuphantomlib.units_util import set_as_quantity
from clmuphantomlib.log import say, is_verbose
from clmuphantomlib.geometry import get_dist2_between_2pt, get_closest_pt_on_line, get_dist2_from_pt_to_line_nb, get_ray_unit_vec, get_rays_unit_vec
from clmuphantomlib.sph_interp import get_sph_interp, get_h_from_rho, get_no_neigh
from clmuphantomlib.units_util import set_as_quantity, set_as_quantity_temperature, get_units_field_name
from clmuphantomlib.eos.base import EoS_Base
#  import (general)
import numpy as np
from numpy import typing as npt
import numba
from numba import jit, prange
import sarracen
import itertools

# settings
#
#   imported from script_input.py file

from script_PhLocAxes__input import interm_dir, verbose, PHOTOSPHERE_TAU, JOB_PROFILES
from _sharedFuncs import mpdf_read


# set metadata
with open("_metadata__input.json", 'r') as f:
    metadata = mupl.json_load(f)
metadata['Title'] = "Getting photosphere size on x, y, z axes"
metadata['Description'] = f"""Tracing 6 rays on +x, -x, +y, -y, +z, -z directon and get photosphere size, h, rho, u, T from them."""

plt.rcParams.update({'font.size': 20})

def get_photosphere_on_ray(
    pts_on_ray            : np.ndarray,
    dtaus                 : np.ndarray,
    pts_order             : np.ndarray,
    sdf                   : sarracen.SarracenDataFrame,
    ray                   : np.ndarray,
    calc_params           : list       = ['loc', 'R1'],
    hfact                 : float      = None,
    mpart                 : float      = None,
    eos                   : EoS_Base   = None,
    sdf_units             : dict       = None,
    ray_unit_vec          : np.ndarray = None,
    kernel                : sarracen.kernels.base_kernel = None,
    do_skip_zero_dtau_pts : bool       = True,
    do_use_precise_calc   : bool       = True,
    photosphere_tau       : float      = 1.,
    waypts_max_tau        : None|float = None, #746.,    # exp(-746) is at float64 machine precision limit
    waypts_nsample        : int        = 1024,
    ndim                  : int        = 3,
    return_as_quantity    : bool|None  = True,
    verbose : int = 3,
) -> tuple[dict, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Calc the location where the photosphere intersect with the ray.

    Assuming 3D.

    
    Parameters
    ----------
    pts_on_ray, dtaus, pts_order
        output from get_optical_depth().

        pts_on_ray: np.ndarray
            Orthogonal projections of the particles' locations onto the ray.
        
        dtaus: np.ndarray
            Optical depth tau contributed by each particles. In order of the original particles order in the dump file.
            Remember tau is a dimensionless quantity.
        
        pts_order: np.ndarray
            indices of the particles where dtaus are non-zero.
            The indices are arranged by distances to the observer, i.e. the particles closest to the observer comes first, 
            and the furtherest comes last.

    sdf: sarracen.SarracenDataFrame
        Must contain columns: x, y, z, h    # kappa, rho,
        
    ray: (2, 3)-shaped numpy array, i.e. [pt1, pt2]
        2 points required to determine a line.
        The line is described as X(t) = pt1 + t*(pt2-pt1)
        First  point pt1 is the reference of the distance if R1 is calc-ed.
        Second point pt2 points to the observer, and is closer to the observer.

    calc_params: list or tuple of str
        parameters to be calculated / interpolated at the photosphere location.
        Results will be put into the photosphere dict in the output.
        Acceptable input: (Note: will always calc 'loc' if calc_params is not empty)
            'is_found': will return bool.
                Will always be outputted regardless of in calc_params or not.
            'loc': will return (3,)-shaped numpy array.
                photophere location.
            'R1' : will return float.
                distance between photosphere location and the ray[0].
                Could be negative if loc is on the other side of the ray.
            'nneigh': will return int.
                Number of neighbour particles of the photosphere loc.
            'rho': will return float.
                density at the photosphere.
            'u': will return float.
                specific internel energy at the photosphere.
            'h'  : will return float.
                smoothing length at the photosphere.
                Will always calc 'rho' if to calc 'h'.
            'T'  : will return float.
                Temperature at the photosphere.
                Warning: if not supplied 'temp' keyword in sdf_units, will return in cgs units.
    
    hfact, mpart: float
        Only useful if you are calc-ing 'h'
        $h_\\mathrm{fact}$ and particle mass used in the phantom sim.
        If None, will get from sdf.params['hfact'] and sdf.params['mass']

    eos: .eos.base.EoS_BASE
        Only useful if you are calc-ing 'T'
        Equation of state object defined in eos/base.py

    sdf_units: dict
        Only useful if you are calc-ing 'T'
        in which case, supply rho, u, and T units in this dict
        e.g.
        sdf_units = {
            'density': units.Msun / units.Rsun**3,
            'specificEnergy': units.Rsun**2 / units.s**2,
            'temp': units.K,
        }
    
    ray_unit_vec: (3,)-shaped np.ndarray
        unit vector of ray. will calc this if not supplied.
        
    kernel: sarracen.kernels.base_kernel
        Smoothing kernel for SPH data interpolation.
        If None, will use the one in sdf.
        
    do_skip_zero_dtau_pts: bool
        Whether or not to skip particles with zero dtaus (i.e. no contribution to opacity) to save computational time.
        If skiped, these particles' locs will be excluded from results as well

    do_use_precise_calc: bool
        If True, will use the new tau interpolation algorithm used in LCGen.
        otherwise, will just interpolate from dtaus input (i.e. assuming each SPH particle is a point mass)
        
    photosphere_tau: float
        At what optical depth (tau) is the photosphere defined.

    waypts_max_tau: float
        maximum optical depth that we will look inside.

    return_as_quantity: bool | None
        If True or None, the results in photosphere will be returned as a astropy.units.Quantity according to sdf_units.
        (pts_waypts, pts_waypts_t, taus_waypts) will also be returned as numpy array and NOT as Quantity.
        The diff between True and None is that True will raise an error if units not supplied in sdf_units,
        while None will just return as numpy array in such case.
        
    
    verbose: int
        How much warnings, notes, and debug info to be print on screen. 


    Returns
    -------
    photosphere, (pts_waypts, pts_waypts_t, taus_waypts)

    photosphere: dict
        dict of values found at the photosphere intersection point with the ray.
        will always have 

    pts_waypts: (npart, 3)-shaped numpy array
        location of the waypoints on ray

    pts_waypts_t: (npart)-shaped numpy array
        distance of the waypoints from ray[0]

    taus_waypts: (npart)-shaped numpy array
        optical depth at the waypoints.

    kcs_waypts: (npart)-shaped numpy array
        summed column kernel value at the waypoints.
        
    """

    # init
    ray = np.array(ray)
    if ray_unit_vec is None:
        ray_unit_vec = get_ray_unit_vec(ray)
    if kernel is None:
        kernel = sdf.kernel
    kernel_col, kernel_csz, _, _ = get_col_kernel_funcs(kernel)
    if do_skip_zero_dtau_pts:
        pts_order = pts_order[np.where(dtaus[pts_order])]
    ray_0 = np.asarray(ray[0])
    pts_ordered    = np.asarray(sdf[['x', 'y', 'z']].iloc[pts_order])
    hs_ordered     = np.asarray(sdf[ 'h'           ].iloc[pts_order])
    #kappas_ordered = np.array(sdf[ 'kappa'       ].iloc[pts_order])
    #rhos_ordered   = np.array(sdf[ 'rho'         ].iloc[pts_order])
    dtaus_ordered  = dtaus[pts_order]
    pts_on_ray_ordered = pts_on_ray[pts_order]
    npart_ordered = pts_ordered.shape[0]


    
    # get waypts (way points) for pts (point locations) and taus (optical depths)
    #  waypts are suitable for linear interpolation
    #  taus_waypts[0] is 0; taus_waypts[-1] is total optical depth from the object

    
    #  step 1: determine the waypts location by assuming pts as balls with constant kappa and density
    
    #   step 1a: getting the size of pts balls on the ray
    pts_dist2_to_ray = get_dist2_between_2pt(pts_ordered, pts_on_ray_ordered)
    #    Assuming a h radius ball
    pts_radius = kernel.get_radius() * hs_ordered
    pts_size_on_ray = pts_radius**2 - pts_dist2_to_ray
    # put a small number (1e-8*h) in case of negative pts_size_on_ray, so that the code does not freak out
    pts_size_on_ray_min = 1e-8*hs_ordered
    pts_size_on_ray = np.where(pts_size_on_ray < pts_size_on_ray_min**2, pts_size_on_ray_min, pts_size_on_ray**0.5)
    #pts_size_on_ray = dtaus_ordered / (kappas_ordered * rhos_ordered)    # does not work because rho is not a constant within the particle

    #   step 1b: getting the waypoint locs
    pts_on_ray_t_ordered = np.sum((pts_on_ray_ordered - ray_0) * ray_unit_vec, axis=-1)
    if waypts_max_tau is None:
        # calc waypoints by interpolating ALL relavent particle's positions
        #    pts_waypts_t: the distance from waypts to ray_0 (negative if in the opposite direction)
        pts_waypts_t = np.interp(    # 5 data points in between pts, so npart_ordered*5 - 1
            np.linspace(0, npart_ordered-1, (npart_ordered-1)*5 + 1),
            np.linspace(0, npart_ordered-1, npart_ordered),
            pts_on_ray_t_ordered)
        pts_waypts_t = np.concatenate((    # add before and after the first particles
            np.linspace(pts_on_ray_t_ordered[0] + pts_radius[0], pts_on_ray_t_ordered[0], 10),
            pts_waypts_t,
            np.linspace(pts_on_ray_t_ordered[-1] - pts_radius[-1], pts_on_ray_t_ordered[-1], 10),
        ))
    else:
        # calc waypoints by interpolating between tau = 0 and tau = waypts_max_tau
        ind = min(np.searchsorted(np.cumsum(dtaus_ordered), waypts_max_tau), npart_ordered-1)
        pts_waypts_t = np.linspace(
            pts_on_ray_t_ordered[0  ] + pts_radius[0  ],
            pts_on_ray_t_ordered[ind] - pts_radius[ind],
            waypts_nsample+1, endpoint=False)[1:]    # remove head point to avoid division by zero warnings
    pts_waypts = ray_0 + pts_waypts_t[:, np.newaxis] * ray_unit_vec[np.newaxis, :]
        

    #   step 1c: sort waypoint locs
    #    sorting should not be necessary, but just in case
    pts_waypts_inds = np.argsort(pts_waypts_t)[::-1]
    pts_waypts   = pts_waypts[  pts_waypts_inds]
    pts_waypts_t = pts_waypts_t[pts_waypts_inds]
    
    #  step 2: determine the waypts optical depth
    taus_waypts = np.zeros(pts_waypts.shape[0])
    kcs_waypts  = np.zeros_like(taus_waypts)
    if do_use_precise_calc:
        # re-calc
        mkappa_div_h2_ordered = np.asarray(sdf['m'].iloc[pts_order] * sdf['kappa'].iloc[pts_order] / hs_ordered**2)
        for j in range(npart_ordered):
            h = hs_ordered[j]
            q_xy_j = pts_dist2_to_ray[j]**0.5 / h
            t_j = pts_on_ray_t_ordered[j]
            dkcs_waypts = np.asarray([
                kernel_csz(q_xy_j, -(t_k - t_j)/h, ndim) for t_k in pts_waypts_t])
            kcs_waypts  += dkcs_waypts
            taus_waypts += mkappa_div_h2_ordered[j] * dkcs_waypts
    else:
        # interpolate from given dtau (not the same as in LCGen)
        for waypt_t, h, dtau in zip(pts_waypts_t, hs_ordered, dtaus_ordered):
            hr = h * kernel.get_radius()
            # Note: np.interp assumes xp increasing, so we need to reverse this
            taus_waypts += np.interp(pts_waypts_t[::-1], [waypt_t-hr, waypt_t+hr], [dtau, 0.], left=dtau, right=0.)[::-1]
        

    # prepare answers
    # is found?
    if not taus_waypts.size:
        taus_max = 0
    elif np.isfinite(taus_waypts[-1]):
        # in case there is nan in the later part of the array
        taus_max = taus_waypts[-1]
    else:
        taus_max = np.nanmax(taus_waypts)
    photosphere = {
        'is_found': (taus_max > photosphere_tau)
    }
    
    # get photosphere parameters
    if calc_params:
        # always calc location if anything needs to be calc-ed
        photosphere['loc'] = np.array([
            np.interp(photosphere_tau, taus_waypts, pts_waypts[:, ax], right=np.nan) if taus_waypts.size else np.nan
            for ax in range(pts_waypts.shape[1])
        ])

        # do prerequisite check
        calc_params = list(calc_params)
        if 'h' in calc_params:
            if 'rho' not in calc_params: calc_params.append('rho')
        if 'T'   in calc_params:
            if 'rho' not in calc_params: calc_params.append('rho')
            if 'u'   not in calc_params: calc_params.append('u')

        # first calc prerequisites
        calc_these = []
        for calc_name in calc_params:
            if   calc_name == 'loc':
                # already calc-ed
                pass
            elif calc_name == 'R1':
                photosphere['R1']  = np.interp(photosphere_tau, taus_waypts, pts_waypts_t, right=np.nan) if taus_waypts.size else np.nan
            elif calc_name in {'rho', 'u'}:
                photosphere[calc_name]  = get_sph_interp(sdf, calc_name, photosphere['loc'], kernel=kernel, verbose=verbose)
            elif calc_name in {'nneigh'}:
                photosphere[calc_name]  = get_no_neigh(sdf, photosphere['loc'], kernel=kernel, verbose=verbose)
            else:
                calc_these.append(calc_name)
    
        # now the rest
        for calc_name in calc_these:
            if calc_name == 'h':
                if hfact is None: hfact = sdf.params['hfact']
                if mpart is None: mpart = sdf.params['mass']
                photosphere['h']  = get_h_from_rho(photosphere['rho'], mpart, hfact)
                if np.isnan(photosphere['h']):
                    # fix for APR case where mpart is nan
                    photosphere['h'] = get_sph_interp(sdf, 'h', photosphere['loc'], kernel=kernel, verbose=verbose)
            elif calc_name == 'T':
                if eos   is None: raise ValueError("get_photosphere_on_ray(): Please supply equation of state to calculate temperature.")
                try:
                    photosphere['T']  = eos.get_temp(
                        set_as_quantity(photosphere['rho'], sdf_units['density']),
                        set_as_quantity(photosphere['u'  ], sdf_units['specificEnergy']))
                    if 'temp' in sdf_units:
                        photosphere['T'] = set_as_quantity_temperature(photosphere['T'], sdf_units['temp']).value
                    else:
                        photosphere['T'] = photosphere['T'].value
                except ValueError:
                    # eos interp could go out of bounds if it's a tabulated EoS
                    # which will raise a Value Error
                    photosphere['T'] = np.nan
            else:
                # just interpolate it (#IT-JUST-WORKS)
                photosphere[calc_name]  = get_sph_interp(sdf, calc_name, photosphere['loc'], kernel=kernel, verbose=verbose)

        # add units
        if return_as_quantity or return_as_quantity is None:
            for calc_name in photosphere.keys():
                if calc_name not in {'is_found', 'nneigh'}:
                    # find appropriate unit
                    try:
                        unit_field_name = get_units_field_name(calc_name)
                    except NotImplementedError:
                        # failed to find unit type
                        unit_field_name = None
                    if unit_field_name in sdf_units.keys():
                        # add units
                        photosphere[calc_name] = set_as_quantity(photosphere[calc_name], sdf_units[unit_field_name])
                    # errors
                    elif unit_field_name is None:
                        if is_verbose(verbose, 'warn'):
                            say('warn', None, verbose,
                                f"Cannot find the corresponding unit for {calc_name}. Will return as numpy array instead.")
                    elif return_as_quantity is not None:
                        raise ValueError(f"Please supply {unit_field_name} in sdf_units.")
        
        
    return photosphere, (pts_waypts, pts_waypts_t, taus_waypts, kcs_waypts)



def write_ph_loc_axes(
    job_profile : dict,
    file_indexes : np.ndarray,
    rays_dir_def : dict,    # dict of list
    eoses : (mupl.eos.base.EoS_Base, mupl.eos.mesa.EoS_MESA_opacity),
    calc_params = ['loc', 'R1', 'rho', 'h', 'T', 'kappa', 'kappaDust', 'srcfunc', 'Skapparho'],
    photosphere_tau = PHOTOSPHERE_TAU,
    verbose : int = 2,
    interp_method: str = 'basic',    # 'basic' or 'improved'
):

    """Writing the photosphere locations of each dump to json files.

    Notes:
    Using mpdf.params['hfact']
    """
    
    #mpdf = mupl.MyPhantomDataFrames()

    
    job_name = job_profile['job_name']
    #X = job_profile['X']
    #ieos = job_profile['ieos']

    eos, eos_opacity = eoses

    
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
            'mpart_Msun' : None,
            'data': {},
            'rays_dir': rays_dir_def,
            'rays': {},
        }  
        for key in rays_dir.keys():
            photosphere_pars['data'][key] = {}

        # read data
        mpdf = mpdf_read(job_name, file_index, eos_opacity, mpdf=None, reset_xyz_by='CoM', do_extrap=False, verbose=verbose)
        sdf = mpdf.data['gas']
        sdf['srcfunc'] = mpdf.const['sigma_sb'] * sdf['T']**4 / pi
        sdf['Skapparho'] = sdf['srcfunc'] * sdf['kappa'] * sdf['rho']
        
        #mpdf.read(job_name, file_index, reset_xyz_by='CoM', verbose=verbose)
        #if 'Tdust' in mpdf.data['gas'].columns:
        #    mpdf.data['gas']['T'] = mpdf.data['gas']['Tdust']
        #elif 'temperature' in mpdf.data['gas'].columns:
        #    mpdf.data['gas']['T'] = mpdf.data['gas']['temperature']
        #if 'kappa' not in mpdf.data['gas'].keys():
        #    # get kappa from mesa table in cgs units
        #    mpdf.data['gas']['kappa'] = eos_opacity.get_kappa(
        #        mpdf.get_val('rho', copy=False),
        #        mpdf.get_val('T', copy=False),
        #        do_extrap=True,
        #        return_as_quantity=False)
        ## translate to phantom units
        #mpdf.calc_sdf_params(
        #    calc_params=['kappa',], #'R1',
        #    calc_params_params={'ieos': ieos, 'X':X, 'overwrite':False, 'kappa_translate_from_cgs_units':True},
        #    verbose=verbose,
        #)
        hfact = mpdf.params['hfact']
        mpart = mpdf.params['mass'] if 'apr_level' not in sdf else np.nan
        
        photosphere_pars['time_yr'] = mpdf.get_time().to_value(units.year)
        photosphere_pars['orbsep_Rsun'] = mpdf.get_orb_sep().to_value(units.Rsun)
        photosphere_pars['mpart_Msun']  = (mpart * mpdf.units['mass']).to_value(units.Msun)


        ## construct kdtree (since we are not changing x, y, z label here)
        #sdf_all_kdtree = kdtree.KDTree(np.array(sdf[['x', 'y', 'z']], copy=False))

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
            if verbose:
                debug_info(
                    'write_ph_loc_axes()', verbose,
                    f"{np.count_nonzero(sdf_selected_indices)} particles are close enough to the ray to have effects."
                )
            sdf = sdf.iloc[sdf_selected_indices]
            pts = np.array(sdf[['x', 'y', 'z']])    # (npart, 3)-shaped array


            # get optical depth
            if verbose:
                debug_info(
                    'write_ph_loc_axes()', verbose,
                    f"{ray = }"
                )
            pts_on_ray, dtaus, pts_order = mupl.light.get_optical_depth_by_ray_tracing_3D(sdf=sdf, ray=ray)
            photosphere, (pts_waypts, pts_waypts_t, taus_waypts, kcs_waypts) = get_photosphere_on_ray(
                pts_on_ray, dtaus, pts_order, sdf, ray,    # remove dtaus to force recalc in LCGen way
                calc_params = calc_params,
                hfact = hfact, mpart=mpart, eos=eos, sdf_units=mpdf.units,
                ray_unit_vec=ray_unit_vec, verbose=verbose, photosphere_tau=photosphere_tau,
                return_as_quantity=False,
            )

            if is_verbose(verbose, 'debug'):
                say('debug', None, verbose,
                    f"{photosphere = }\n",
                    f"{len(taus_waypts) = }\n",
                )
            
            photosphere_pars['data'][key] = photosphere
            photosphere_pars['data'][key]['size'] = photosphere['R1']
            ph_d = photosphere_pars['data'][key]

            # save interpolated data on ray at waypoints
            ph_d[ 'R1_on_ray'] = pts_waypts_t
            ph_d['loc_on_ray'] = ray[0][np.newaxis, :] + ph_d['R1_on_ray'][:, np.newaxis] * ray_unit_vec[np.newaxis, :]
            ph_d['tau_on_ray'] = np.interp(ph_d['R1_on_ray'], pts_waypts_t[::-1], taus_waypts[::-1])
            ph_d[ 'kc_on_ray'] = np.interp(ph_d['R1_on_ray'], pts_waypts_t[::-1],  kcs_waypts[::-1])
            ph_d['rho_on_ray'] = mupl.sph_interp.get_sph_interp(sdf, 'rho', ph_d['loc_on_ray'], verbose=verbose, method=interp_method)
            ph_d[  'h_on_ray'] = mupl.sph_interp.get_h_from_rho(ph_d['rho_on_ray'], mpart=mpart, hfact=hfact)
            # ph_d[  'u_on_ray'] = mupl.sph_interp.get_sph_interp(
            #     sdf, 'u'  , ph_d['loc_on_ray'], verbose=verbose, method=interp_method)
            # ph_d[  'T_on_ray'] = eos.get_temp(
            #     set_as_quantity(photosphere_pars['data'][key]['rho_on_ray'], mpdf.units['density']),
            #     set_as_quantity(photosphere_pars['data'][key]['u_on_ray']  , mpdf.units['specificEnergy']),
            #     return_as_quantity=False, bounds_error=False)
            for par in calc_params:
                if par in {'R1', 'loc', 'tau', 'rho', 'h',}:    # 'u', 'T',
                    continue    # already handled
                elif par in sdf:    # just interpolate it    #IT-JUST-WORKS
                    ph_d[f'{par}_on_ray'] = mupl.sph_interp.get_sph_interp(
                        sdf, par, ph_d['loc_on_ray'], verbose=verbose, method=interp_method)
                else:
                    say('warn', None, verbose, f"{par} not found in data dump '{job_name}_{file_index:05}'.")
            photosphere_pars['data'][key]['nneigh_on_ray']=mupl.sph_interp.get_no_neigh(
                sdf, ph_d['loc_on_ray'], kernel_rad=kernel_radius)

            # save actual relevant particle data (pts)
            #    Rt_pts: equivalant of R1_on_ray but for particles
            ph_d[ 'Rt_at_pts'] = np.sum((pts_on_ray[pts_order] - np.asarray(ray[0])) * ray_unit_vec, axis=-1)
            ph_d['Rxy_at_pts'] = get_dist2_between_2pt(pts[pts_order], pts_on_ray[pts_order])**0.5
            ph_d['tau_at_pts'] = np.cumsum(dtaus[pts_order])
            for par in calc_params:
                if par in {'R1', 'loc', 'tau',}:
                    continue
                elif par in sdf:
                    ph_d[f'{par}_at_pts'] = np.asarray(sdf[par].iloc[pts_order])
                else:
                    say('warn', None, verbose, f"{par} not found in data dump '{job_name}_{file_index:05}'.")


        hdf5_dump(
            photosphere_pars,
            f"{interm_dir}{job_profile['nickname']}_{file_index:05d}.photospherePars.xyz.hdf5",
            metadata=metadata)
        del mpdf

    return photosphere_pars
