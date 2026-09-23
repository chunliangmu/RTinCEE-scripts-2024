# Scripts for Mu et al., 2026

Scripts for the paper *[Dust Formation in Common Envelope Binary Interactions — III. Lightcurves](https://doi.org/10.48550/arXiv.2606.26495)*.

> Author: Chunliang Mu  
> Requrie python 3.10+
>
> These scripts are written by me for my PhD project **Radiative Transfer (RT) in Common Envelope Evolution (CEE)** (a.k.a. "Non-adiabatic Common Envelope Simulation of Massive Stars").  
> Creator: ***Chunliang Mu*** (PhD student at Macquarie University 2023-2026)  
> Principal Supervisor: Professor Orsola De Marco  
> Associate Supervisor: Professor Mark Wardle  

For the code to work properly, put my `clmuphantomlib` in the directory where the symbolic link is pointing to. (See github link at the bottom of this page.)
My personal set up has this repository's files in folder `[PATH]/scripts/`, and the `clmuphantomlib` repository's files in folder `[PATH]/src/clmuphantomlib/`.

**Note: Please cite the sarracen paper if you use this code (see below link for the sarracen repository description), since this code uses sarracen behind the scene.**

## Dependencies

- Python libraries:
	- `python3` (version >= 3.10)
	- `numpy scipy astropy h5py numba matplotlib ipympl moviepy`
	- [`sarracen`](https://github.com/ttricco/sarracen)
    - [`clmuphantomlib`](https://github.com/chunliangmu/clmuphantomlib)

## Externel files

- `.gitignore`: obtained from https://github.com/github/gitignore/blob/main/Python.gitignore under CC0-1.0 license.

## Useful links

- `phantom` [GitHub](https://github.com/danieljprice/phantom)
- `sarracen` [GitHub](https://github.com/ttricco/sarracen)
- `clmuphantomlib` [GitHub](https://github.com/chunliangmu/clmuphantomlib)
