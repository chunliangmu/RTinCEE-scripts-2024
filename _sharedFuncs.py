#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Shared helper functions
"""

# imports and internal settings
import numpy as np
from numpy import typing as npt
from numpy import pi
from astropy import units
from astropy import constants as const
from clmuphantomlib.log import is_verbose, say
from clmuphantomlib import MyPhantomDataFrames, get_eos, get_eos_opacity
from clmuphantomlib.eos   import EoS_MESA_opacity
from clmuphantomlib.units_util import get_val_in_unit, get_units_field_name



# - constants -


_UNITS_OUT: dict[str, units.Unit] = {
    'time': units.yr,
    'energy': units.erg, #units.Lsun*units.yr,
    'linearMomentum': None,
    'angularMomentum': None,
    'mass': units.Msun,
    'temp': units.K,
}

_DTYPE_LUIS_2MD = [
    (('01        time', 'time'      ), np.float64),
    (('02        ekin', 'E_kin'     ), np.float64),
    (('03      etherm', 'E_hea'     ), np.float64),
    (('04        emag', 'E_mag'     ), np.float64),
    (('05        epot', 'E_pot'     ), np.float64),
    (('06        etot', 'E_tot'     ), np.float64),
    (('07        erad', 'E_rad'     ), np.float64),
    (('08      totmom', 'p_tot'     ), np.float64),
    (('09      angtot', 'l_tot'     ), np.float64),
    (('10     rho max', 'rho_max'   ), np.float64),
    (('11     rho ave', 'rho_ave'   ), np.float64),
    (('12          dt', 'dt'        ), np.float64),
    (('13   totentrop', 'totentrop' ), np.float64),
    (('14     rmsmach', 'rmsmach'   ), np.float64),
    (('15        vrms', 'vrms'      ), np.float64),
    (('16        xcom', 'xcom'      ), np.float64),
    (('17        ycom', 'ycom'      ), np.float64),
    (('18        zcom', 'zcom'      ), np.float64),
    (('19   alpha max', 'alpha_max' ), np.float64),
    (('20    temp max', 'T_max'     ), np.float64),
    (('21    temp ave', 'T_ave'     ), np.float64),
    (('22    temp min', 'T_min'     ), np.float64),
]



# --- mpdf reading ---


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
    use_Tscales: None|str = None,
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
        filename = f"{job_name}_T{use_Tscales}.npy"
        scales = np.load(filename)
        # loop over iorig to make sure indexes mismatch doesn't happen
        # delete negative Tscales
        inds_del = scales['T_scale'] < 0
        mpdf.data['gas'].drop([j-1 for j in scales['iorig'][inds_del]], inplace=True, errors='ignore')
        sdf = mpdf.data['gas']
        failed_count = 0
        for j, T_scale in zip(scales['iorig'][~inds_del], scales['T_scale'][~inds_del]):
            #    note: particles can be deleted, so we need to try
            try:
                assert sdf.at[j-1, 'iorig'] == j
                sdf.at[j-1, 'T'] *= T_scale
            except KeyError:
                failed_count += 1
        say('note', None, verbose,
            f"Using '{filename}' to scale temperature",
            f"out of {len(scales)} particles selected, {np.count_nonzero(inds_del)} are to be ignored in post-processing,"
            f"{failed_count} out of the rest {len(scales) - np.count_nonzero(inds_del)} particles not found.")
        
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
            say('debug', None, verbose,
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
    method : str  = 'scale',
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

    method: {'scale', 'cut', 'delete'}

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

    inds_outer = mpdf.get_val('R1')[inds] - mpdf.get_val('h')[inds] > R_ph
    if method == 'cut':
        scales['T_scale'][inds_outer] = -1.
    elif method == 'delete':
        scales['T_scale'][~inds_outer] = 1.
        scales['T_scale'][inds_outer] = -1.

    scales.sort(order='iorig')

    if do_save:
        filename = f"{job_name}_T{method}.npy"
        say('note', None, verbose, f"Saving to {filename}")
        with open(filename, 'wb') as fp:
            np.save(fp, scales)

    say('note', None, verbose,
        f"Scaling {len(scales)} particles out of {len(mpdf.data['gas'])}",
        f"\t({len(scales) / len(mpdf.data['gas']) * 100:.2f}%)",
        f"{min(scales['T_scale']) = :5.3f}    {max(scales['T_scale']) = :5.3f}",
        f"\tDeleting {np.count_nonzero(scales['T_scale'] == -1.)} particles;",
        f"\tLeaving {np.count_nonzero(scales['T_scale'] == 1.)} particles untouched.",
    )

    return scales



# --- .ev files reading ---


def _get_dtype_code_from_header(
    txt: str,
    sep: str = "]   [",
    num_len: int = 2,
    do_print:bool=True,
) -> str:
    """Translate .ev file header into numpy dtype
    
    e.g. something from
    
    txt="[01        time]   [02        ekin]   [03      etherm]   [04        emag]   [05        epot]   [06        etot]   [07        erad]   [08      totmom]   [09      angtot]"
    
    into

    ans="(('01        time', 'time'   ), np.float64),
    (('02        ekin', 'ekin'   ), np.float64),
    (('03      etherm', 'etherm' ), np.float64),
    (('04        emag', 'emag'   ), np.float64),
    (('05        epot', 'epot'   ), np.float64),
    (('06        etot', 'etot'   ), np.float64),
    (('07        erad', 'erad'   ), np.float64),
    (('08      totmom', 'totmom' ), np.float64),
    (('09      angtot', 'angtot' ), np.float64),"

    num_len is the length of the numbers ('2' -> 1, '02' -> 2, '002' -> 3)
    
    """
    txt_list = txt.split(sep)
    txt_list[ 0] = txt_list[ 0].split("[")[-1]
    txt_list[-1] = txt_list[-1].split("]")[ 0]

    names = ['_'.join(entry[num_len:].strip().split(' ')) for entry in txt_list]
    names_maxlen = np.max([len(name) for name in names])

    ans_list = [

        f"(('{entry}', '{name}' {' '*(names_maxlen-len(name))}), np.float64),"
        for entry, name in zip(txt_list, names)
    ]

    ans = '\n'.join(ans_list)
    if do_print: print(ans)

    return ans



def pa_read_energy(
    filepath,
    units_in : None | dict[str, units.Unit] = None,
    units_out: dict[str, units.Unit] = _UNITS_OUT.copy(),
) -> tuple[npt.NDArray, dict[str, npt.NDArray|units.Quantity]]:
    """Read phantom analysis data - CE - energy.
    
    units_out is only used if units_in is not None

    Example usage:
        _, d = pa_read_energy(filepath_ev, units_in=mpdf.units)
    """
    data = np.genfromtxt(filepath, dtype=[
        ((' 1         time', 'time'),  np.float64),
        ((' 2total energy ', 'E_tot'), np.float64),
        ((' 3  pot energy ', 'E_pot'), np.float64),
        ((' 4  kin energy ', 'E_kin'), np.float64),
        ((' 5therm energy ', 'E_hea'), np.float64),
        ((' 6    sink pot ', 'E_pot_sink'), np.float64),    # "does not include sink-gas potential energy"
        ((' 7    sink kin ', 'E_kin_sink'), np.float64),
        ((' 8    sink orb ', 'E_orb_sink'), np.float64),    # "sink kin + sink pot"
        ((' 9    comp orb ', 'E_orb_comp'), np.float64),
        (('10     env pot ', 'E_pot_env' ), np.float64),
        (('11  env energy ', 'E__env'), np.float64),
        (('12   bound kin ', 'E_kin_bound'), np.float64),
        (('13 unbound kin ', 'E_kin_unbound'), np.float64),
        (('14  bound mass ', 'm_bound'), np.float64),
        (('15unbound mass ', 'm_unbound'), np.float64),
        (('16     p-p pot ', 'E_pot_pp'), np.float64),
        (('17     p-s pot ', 'E_pot_ps'), np.float64),
        (('18 tot ang mom ', 'l_tot'), np.float64),
        (('19   b ang mom ', 'l_bound'), np.float64),
        (('20  ub ang mom ', 'l_unbound'), np.float64),
        (('21 orb ang mom ', 'l_orb'), np.float64),
        (('22  gas energy ', 'E__gas'), np.float64),
        (('23    fallback ', 'm_fallback'), np.float64),
        (('24fallback mom ', 'l_fallback'), np.float64),
    ])

    # the following has the same units per entry
    grp_by_units: dict[str, npt.NDArray|units.Quantity] = {
        't': data[['time']],
        'E': data[[    # having the same energy unit
            'E_tot', 'E_pot', 'E_kin', 'E_hea',
            'E_pot_sink', 'E_kin_sink', 'E_orb_sink', 'E_orb_comp', 'E_pot_env',
            'E__env', 'E_kin_bound', 'E_kin_unbound','E_pot_pp', 'E_pot_ps', 'E__gas']],
        'l': data[[    # angularMomentum
            'l_tot', 'l_bound', 'l_unbound', 'l_orb', 'l_fallback']],
        'm': data[[    # mass
            'm_bound', 'm_unbound', 'm_fallback']],
    }

    for label in grp_by_units.keys():
        field = get_units_field_name(label)
        if units_in  is not None and field in units_in  and units_in[field]  is not None:
            grp_by_units[label] *= units_in[field]
        if units_out is not None and field in units_out and units_out[field] is not None:
            grp_by_units[label] = grp_by_units[label].to(units_out[field])
    
    return data, grp_by_units



def pa_read_ev(
    filepath,
    dtype : npt.DTypeLike = _DTYPE_LUIS_2MD,
    grp_by_units_dict: None | dict[str, list[str]] = {
        't': ['time'],
        'E': [    # having the same energy unit
            'E_kin', 'E_hea',  'E_mag', 'E_pot', 'E_tot', 'E_rad',],
        'p': [    # Momentum
            'p_tot'],
        'l': [    # momentum
            'l_tot',],
        'T': [    # Temperature
            'T_max', 'T_ave', 'T_min',],
    },
    units_in : None | dict[str, units.Unit] = None,
    units_out: dict[str, units.Unit] = _UNITS_OUT.copy(),
) -> tuple[npt.NDArray, dict[str, npt.NDArray|units.Quantity]]:
    """Read phantom ev.
    
    units_out is only used if units_in is not None

    Example usage:
        _, d = pa_read_energy(filepath_ev, units_in=mpdf.units)
    """
    data = np.genfromtxt(filepath, dtype=dtype)

    # the following has the same units per entry
    grp_by_units: dict[str, npt.NDArray|units.Quantity] = {
        k: data[v] for k, v in grp_by_units_dict.items()
    }

    for label in grp_by_units.keys():
        field = get_units_field_name(label)
        if units_in  is not None and field in units_in  and units_in[field]  is not None:
            grp_by_units[label] *= units_in[field]
            if units_out is not None and field in units_out and units_out[field] is not None:
                grp_by_units[label] = grp_by_units[label].to(units_out[field])
    
    return data, grp_by_units
