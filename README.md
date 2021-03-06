# X-regression

**X-regression** library should serve to regression of particular phenomena, e.g. 11-year solar cycle, and their climatic effects. Your feedback is highly appreciated.

## Required package installation
`pip install -r requirements.txt`

## Usage
Generate QBO orthogonal regressors via `eofs_univ_ccmi_xarray.py` script, e.g.:

`python eofs_univ_ccmi_xarray.py -v jra55 u 1960 1960 2009 jra55_uwnd_1960_2009_zm.nc`

When all regressors available, attribution by MLR to particular phehomena can be done  via `lin_reg_univ_ccmi_xarray.py` script, e.g.:

`python lin_reg_univ_ccmi_xarray.py -v jra55 t 1960 1979 2005 jra55_tmp_1960_2009_zm.nc all_trend`

## Possible acknowledgement via Zenodo:
[![DOI](https://zenodo.org/badge/60927195.svg)](https://zenodo.org/badge/latestdoi/60927195)

## References
[Kuchar, A., Ball, W. T., Rozanov, E. V., Stenke, A., Revell, L., Miksovsky, J., Pisoft P., Peter, T. (2017). On the aliasing of the solar cycle in the lower stratospheric tropical temperature. Journal of Geophysical Research: Atmospheres. http://doi.org/10.1002/2017JD026948](http://onlinelibrary.wiley.com/doi/10.1002/2017JD026948/abstract)



