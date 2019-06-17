# Scripts intended for SOLARIS-HEPPA regression experiments
## Required package installation (for Unix and python3 users)
`python3 -m pip install -r requirements.txt`

To install `scipy` required for `statsmodels` you need to install [BLAS](https://en.wikipedia.org/wiki/Basic_Linear_Algebra_Subprograms) and [LAPACK](https://en.wikipedia.org/wiki/LAPACK) (see [here](https://stackoverflow.com/questions/33368261/what-is-the-easiest-way-to-install-blas-and-lapack-for-scipy/33369271)).

## Usage
`python3 solaris_heppa_analysis-script.py config_SOCOL.yaml`
