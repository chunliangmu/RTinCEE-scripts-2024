#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Shared helper functions
"""

# imports and internal settings
import numpy as np
from numpy import pi
from astropy import units
from astropy import constants as const
from clmuphantomlib.log import is_verbose, say
from clmuphantomlib import MyPhantomDataFrames, get_eos, get_eos_opacity
from clmuphantomlib.eos   import EoS_MESA_opacity
from clmuphantomlib.units_util import get_val_in_unit


def mpdf_read(
    job_name: str,
    file_index: int,
    eos_opacity: None|EoS_MESA_opacity = None,
    mpdf: None|MyPhantomDataFrames = None,
    params    : dict = None,
    reset_xyz_by: str='CoM',
    calc_params :list=['vr', 'R1'],
    kappa_gas : units.Quantity = 2e-4*(units.cm**2/units.g),
    kappa_tol : units.Quantity = 1e-7*(units.cm**2/units.g),
    T_cond_oxy: units.Quantity = 1450 * units.K,
    do_extrap: bool = False,
    use_Tscales: bool = False,
    verbose: int = 3,
) -> MyPhantomDataFrames:
    """Read the dump files and get T and kappa,
    
    using a blend of mesa opacity (for high T) and stored kappa value in the dump (for low T; in my project, this col is nucleation dust opacity).
    If no kappa is found, will use default 2e-4 cm2/g for the low T part instead. You can disable this by setting T_cond_oxy to 0K.

    Supply either eos_opacity or params.
    """
    if mpdf is None:
        mpdf = MyPhantomDataFrames()

    mpdf.read(job_name, file_index, reset_xyz_by=reset_xyz_by, calc_params=calc_params, verbose=verbose)

    if eos_opacity is None:
        eos_opacity = get_eos_opacity(ieos=mpdf.ieos, params=params)
    
    temp_key = {'T', 'temperature', 'Tdust'}.intersection(mpdf.data['gas'].keys()).pop()
    mpdf.data['gas']['T'    ] = mpdf.data['gas'][temp_key]
    if use_Tscales:
        # apply the temperature scales we obtained earilier
        filename = f"{job_name}_Tscales.npy"
        scales = np.load(filename)
        inds = mpdf.data['gas']['iorig'].isin(scales['iorig'])    # warning: may cause indexes mismatch
        # make sure that doesn't happen
        assert np.all(mpdf.data['gas'][inds]['iorig'] == scales['iorig'])
        mpdf.data['gas'].loc[inds, 'T'] *= scales['T_scale']
        
    if 'kappa' in mpdf.data['gas'].keys():
        unit_opacity = units.cm**2/units.g    # opacity units in original phantom dumpfiles
        kappa_mesa = eos_opacity.get_kappa(mpdf.get_val('rho'), mpdf.get_val('T'), do_extrap=do_extrap)
        mpdf.data['gas']['kappa_dust'] = mpdf.data['gas']['kappa'] - kappa_gas.to_value(unit_opacity)
        # fix negative opacities
        mpdf.data['gas']['kappa_dust'] = np.where(
            mpdf.data['gas']['kappa_dust'] < kappa_tol.to_value(unit_opacity),
            0.,
            mpdf.data['gas']['kappa_dust'],
        )
        mpdf.data['gas']['kappa'] = np.where(
            mpdf.data['gas']['T'] < T_cond_oxy,
            mpdf.data['gas']['kappa'],
            mpdf.data['gas']['kappa_dust'] + kappa_mesa.to_value(unit_opacity),
        )
        mpdf.data['gas']['kappa'] = get_val_in_unit(mpdf.data['gas']['kappa'], unit_opacity, mpdf.units['opacity'])
        mpdf.data['gas']['kappa_dust'] = get_val_in_unit(mpdf.data['gas']['kappa_dust'], unit_opacity, mpdf.units['opacity'])
        
        if is_verbose(verbose, 'debug'):
            say('debug', 'mpdf_read()', verbose,
                f"{np.count_nonzero(mpdf.data['gas']['kappa_dust']) = }",
                f"{np.count_nonzero(mpdf.data['gas']['kappa'])      = }")
    else:
        #raise NotImplementedError("non-dusty sims (no kappa column in dump files) not yet implemented")

        # just use mesa opacity
        
        kappa_mesa = eos_opacity.get_kappa(mpdf.get_val('rho'), mpdf.get_val('T'), do_extrap=do_extrap)
        mpdf.data['gas']['kappa'] = np.where(
            mpdf.data['gas']['T'] < T_cond_oxy,
            kappa_gas.to_value( mpdf.units['opacity']),
            kappa_mesa.to_value(mpdf.units['opacity']),
        )

    return mpdf



def gen_Tscales(
    job_name: str,
    T_ph : units.Quantity,
    R_ph : units.Quantity,
    do_save: bool = True,
    params : dict = {'X':0.7, 'Z':0.0},
    verbose: int = 3,
) -> MyPhantomDataFrames:
    """Generate scales at t=0 to scale down temperatures of outter particles.
    
    So temperatures outside the MESA-calculated photosphere
    is scaled down to the MESA values.
    This can be applied to later times,
    so we can see how much the effect of the instable init
    of MESA -> Phantom dump can be.

    ---------------------------------------------------------------------------
    """
    mpdf = mpdf_read(
        job_name, file_index=0, reset_xyz_by='', params=params,
        calc_params=['R1'], verbose=verbose)

    # particle indexes
    inds = mpdf.get_val('R1') + 2*mpdf.get_val('h') > R_ph
    scales = np.zeros(
        np.count_nonzero(inds),
        dtype=[('iorig', np.int64), ('T_scale', np.float64)])
    # particle ids
    #    retrieve back inds from ids by inds_new = mpdf.data['gas']['iorig'].isin(ids)
    scales['iorig'] = mpdf.data['gas']['iorig'][inds]
    scales['T_scale'] = (T_ph / mpdf.get_val('T')[inds]).to_value(
        units.dimensionless_unscaled)
    # scale down only, do not heat the particles
    scales['T_scale'][scales['T_scale'] > 1.] = 1.

    if do_save:
        filename = f"{job_name}_Tscales.npy"
        say('note', None, verbose, f"Saving to {filename}")
        with open(filename, 'wb') as fp:
            np.save(fp, scales)

    return scales
