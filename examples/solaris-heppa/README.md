# Scripts intended for SOLARIS-HEPPA regression experiments
## Required package installation using pip (for Unix and python3.6 users)
`python3 -m pip install -r requirements.txt`

To install `scipy` required for `statsmodels` you need to install [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) (see [here](https://stackoverflow.com/questions/33368261/what-is-the-easiest-way-to-install-blas-and-lapack-for-scipy/33369271)).

## Required package installation using conda
The easiest way to get everything installed is to use [conda](https://conda.io/en/latest/). 

`conda env create --file py36-env.txt`

`source activate py36`

## Usage
### To generate configuration yaml file and run the regression script:
`bash generate_config_refC2.sh 1960 2010 zmta,zmua refC2` runs for for the historical (1960-2010) period, zonal-mean temperature and zonal wind from refC2 simulations
### To run solely the regression script with pregenerated configuration file:
`python3 solaris_heppa_analysis-script.py config_SOCOL3r1i1p1_zmta,zmua_1960-2010.yaml` runs for for the historical (1960-2010) period, zonal-mean temperature and zonal wind from SOCOL3 refC2 simulation
