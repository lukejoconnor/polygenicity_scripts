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
- `estimate_polygenicity.m` reproduces Figure 2 of the paper.

The function which estimates polygenicity from FMR estimates is provided in `FMR` repository. This function in its entirety is:

```matlab
function Pi = compute_polygenicity(x,w,f,finv)
%compute_polygenicity evaluates the function
% Pi_f(x,w) = h^2 / f^-1 (1/h^2 w_1*f(x_1)+...+w_k*f(x_k))
% where h^2 = w_1+...+w_k and f:(0,infty)->(0,infty) is a continuous
% function with inverse finv.
%   For example:
%   \Pi_{entropy}: f = @log and finv=@exp
%   \Pi_{effective}: f=@(x)x and finv=f
%   \Pi_{softmax}: f=@(x)exp(-1./x) and finv=@(x)-1./log(x)
%   \Pi_{softmax} alternative implementation that avoids overflow issues:
%       f=@(x)exp(1/max(x(:)) - 1./x) and finv=@(y,xmax)-1./(-1/xmax + log(y))
%   \Pi_0 [not recommended]: f=@inv and finv=@inv
% 
% Can be applied to matrix-valued x and w (e.g., the jackknife output of FMR),
% in which case it computes polygenicity for each row of the matrix.

h2 = sum(w,2);
w = w./h2;
y = sum(w .* f(x),2);
if nargin(finv) == 1
    Pi = h2 ./ finv(y);
else
    % To avoid overflow issues
    Pi = h2 ./ finv(y, max(x(:)));
end
end
```

## Links and citations
- O'Connor & Sella preprint: https://www.biorxiv.org/content/10.1101/2025.07.10.664154v1
- O'Connor & Sella citation: O'Connor, Luke J., and Guy Sella. "Principled measures and estimates of trait polygenicity." bioRxiv (2025): 2025-07.
- FMR repository: https://github.com/lukejoconnor/FMR
- O'Connor 2021 non-paywalled link: https://www.nature.com/articles/s41588-021-00901-3.epdf?sharing_token=flx8PE5EGIKA7RpaKhSJONRgN0jAjWel9jnR3ZoTv0N7Pnc_k9O2zeCsUKCBmAYoz9yEJzbMB_QL1FfWuvG0UnX1ad9wUpjHqk7ovqIGZcqhYfTjFoKUhoZQNVimIQgn_ZCbpD4IJx18LwQY5QULuXJ6XGkCY30-v-snvrVMFpY%3D
- O'Connor 2021 citation: O’Connor, Luke J. "The distribution of common-variant effect sizes." Nature Genetics 53.8 (2021): 1243-1249.
