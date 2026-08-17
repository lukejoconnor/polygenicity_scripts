# Scripts to reproduce results from O'Connor & Sella 2025

This repository contains scripts to reproduce results from [O'Connor & Sella 2025](https://www.biorxiv.org/content/10.1101/2025.07.10.664154v1). It relies upon the [Fourier Mixture Regression](https://github.com/lukejoconnor/FMR) repository.

## Installation and data download

```bash
git clone https://github.com/lukejoconnor/FMR.git
git clone https://github.com/lukejoconnor/polygenicity_scripts.git
mkdir FMR/matfiles
curl -L -o FMR/matfiles.zip https://www.dropbox.com/sh/mclm1urkxs8ga80/AADDDABQYeGtyQxmom2raMkva
unzip FMR/matfiles.zip -d FMR/matfiles
```

Downloading the data from Dropbox will take several minutes. The size of the download is 6GB. The larger files with LD matrices and Fourier scores are not actually required for these scripts, so you may also choose to download individual files from the [Dropbox link](https://www.dropbox.com/sh/mclm1urkxs8ga80/AADDDABQYeGtyQxmom2raMkva?dl=0).

## Contents

The MATLAB folder contains two scripts:
- `compare_polygenicities.m` reproduces Figure 1 of the paper.
- `estimate_polygenicity.m` reproduces Figure 3 of the paper.

A subdirectory MATLAB/as-is contains a simulation script which was used to produce Figure 2 of the paper. This script is provided as-is; reproducing results yourself will require multiple steps of installation and manual path manipulation.

The Figure 3 analysis uses `MATLAB/helpers/compute_polygenicity.m` to evaluate
Equation 15 directly. Its primary interface is:

```matlab
Pi = compute_polygenicity(sigma2, omega, h2, measure)
```

Here `sigma2` contains the unnormalized FMR component variances in phenotypic-
variance units, each row of `omega` contains component heritability fractions
that sum to one, and `measure` is `entropy`, `effective`, or `softmax`. The
helper retains backward compatibility with the previous `(x,w,f,finv)`
interface. See the function documentation and `analysis/README.md` for the
equation, numerical-stability details, and verification commands.

## Links and citations
- O'Connor & Sella preprint: https://www.biorxiv.org/content/10.1101/2025.07.10.664154v1
- O'Connor & Sella citation: O'Connor, Luke J., and Guy Sella. "Principled measures and estimates of trait polygenicity." bioRxiv (2025): 2025-07.
- FMR repository: https://github.com/lukejoconnor/FMR
- O'Connor 2021 non-paywalled link: https://www.nature.com/articles/s41588-021-00901-3.epdf?sharing_token=flx8PE5EGIKA7RpaKhSJONRgN0jAjWel9jnR3ZoTv0N7Pnc_k9O2zeCsUKCBmAYoz9yEJzbMB_QL1FfWuvG0UnX1ad9wUpjHqk7ovqIGZcqhYfTjFoKUhoZQNVimIQgn_ZCbpD4IJx18LwQY5QULuXJ6XGkCY30-v-snvrVMFpY%3D
- O'Connor 2021 citation: O’Connor, Luke J. "The distribution of common-variant effect sizes." Nature Genetics 53.8 (2021): 1243-1249.
