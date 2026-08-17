# Equation 15 verification

`verify_eq15.py` independently reproduces the published MATLAB calculation
and compares it with the literal Equation 15 implementation. It requires
Python, NumPy, and SciPy; `--published-table` additionally requires openpyxl.

```bash
python analysis/verify_eq15.py \
  --matfiles ../FMR/matfiles \
  --published-table /path/to/NIHMS2174363-supplement-2.xlsx \
  --output analysis/eq15_comparison.csv
```

The required input files are `FMRestimates_32traits.mat`, `FMR_lipids.mat`,
and `name_abbreviations.mat`. The fourth file used by the full Figure 3 MATLAB
script, `32wellpowered_sumstats.mat`, is not needed for the polygenicity
comparison itself.

The verification distinguishes two changes to softmax polygenicity:

1. replacing the normalized-component approximation with literal Equation 15;
2. removing the legacy `1e-256` floor through stable log-sum-exp evaluation.

The latter matters for general risk tolerance and insomnia, whose published
softmax estimates hit the artificial value `-log(1e-256)`.

The verifier fails if the entropy/effective calculations change beyond
floating-point tolerance, if the expected ordering of the three measures is
violated, or (when `--published-table` is supplied) if the legacy calculation
does not reproduce Table S2. `test_compute_polygenicity.m` provides a small
MATLAB fixture for the Equation 15 helper and its backward-compatible API.
