# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + {"toc": true, "cell_type": "markdown"}
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>
# -

import yaml
from pathlib import Path
import argparse
import xarray as xr
import statsmodels.api as sm
import statsmodels.stats.stattools as sms
import numpy as np
import sys
import pandas as pd
import datetime
import dask
from dask.diagnostics import ProgressBar
import dask.multiprocessing
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
#dask.config.set(scheduler='processes')
#dask.config.set(scheduler='threads') 

# +
def _fit_ufunc(y, X, rho = 2):   
    X = X.T
    #print(X.shape)
    #sys.exit()
    nr = X.shape[1]-1
    try:
        mod = sm.GLSAR(y, X, rho, missing = 'drop') # MLR analysis with AR2 modeling
        res = mod.iterative_fit()
        out = np.array((res.params[1:], res.pvalues[1:], \
                        [sms.durbin_watson(res.wresid)]*nr, 
                        [res.rsquared]*nr))
    except:
        out = np.full((4,nr), np.nan)
        
    return out

def deseasonalize(da, var, sel_period = dict(time = slice('1960', '1990'))):
    climatology = da.sel(**sel_period).groupby('time.month').mean('time')
    if var in ['vmro3', 'o3', 'zmo3']:
        anomalies = (da.groupby('time.month') - climatology).groupby('time.month')/climatology*100
    else:
        anomalies = da.groupby('time.month') - climatology
    return anomalies, climatology

dictfilt = lambda x, y: dict([ (i,x[i]) for i in x if i in set(y) ])

def process_out(ds, names, units, config):
    ds['reg'] = names
    ds['var'] = ['coefs', 'p_values', 'DWT', 'CoD']
    ds = ds.to_dataset(dim = 'var')
    
    ds['reg'].attrs['long_name'] = "Regressors's names"
    ds['CoD'].attrs['long_name'] = 'Coefficient of determination'
    ds['DWT'].attrs['long_name'] = 'Durbin-Watson test'
    ds['coefs'].attrs['long_name'] = 'Regression coefficients'
    ds['coefs'].attrs['units'] = units
    ds['p_values'].attrs['long_name'] = 'Statistical significance (p-value)'
    
    now = datetime.datetime.utcnow()
    ds.attrs['description'] = 'Created ' + now.strftime("%Y-%m-%dT%H:%M:%SZ")    
    ds.attrs['config'] = yaml.dump(dictfilt(config, ('analysis',config['analysis']['reg_config'])))
    
    return ds
    
def process_in(da, var, sel_time_dict):
    da = da.sel(**sel_time_dict)
    anomalies, _  = deseasonalize(da, var)
    return anomalies, anomalies['time'], anomalies['month']

def load_regressors(config, sel_time_dict):
    reg_ls = []
    reg_name_ls = []
    reg_dict = config[config['analysis']['reg_config']]
    
    for reg_name in reg_dict.keys():
        reg_name_ls.append(reg_name)
        ds = xr.open_dataset(reg_dict[reg_name]['filename'])[reg_name]
        ds = normalize(ds, reg_dict[reg_name]['norm_type'],  reg_dict[reg_name]['norm_const'])
        reg_ls.append(ds)
        #print(reg_name)

    reg_all = xr.merge(reg_ls).sel(**sel_time_dict)
    
    reg_all['ones'] = (('time'), np.ones(reg_all['time'].shape[0]))
    reg_all = reg_all.to_array()#.T
    reg_all = reg_all.sel(variable = reg_all.coords['variable'].values[::-1])
    
    return reg_all, reg_name_ls

def normalize(da, norm_type, norm_const = None):
    name = da.name

    if norm_type == 'standardize':
        _lambda_ufunc = lambda x: (x - x.mean()) / x.std()
    elif norm_type == 'remove_mean':
        _lambda_ufunc = lambda x: x - x.mean()
    elif norm_type == 'std_scaling':
        _lambda_ufunc = lambda x: x / x.std()
    elif norm_type == 'mean_scaling':
        _lambda_ufunc = lambda x: x / x.mean()
    elif norm_type == 'const_scaling':
        _lambda_ufunc = lambda x: x / norm_const
    elif norm_type == None:
        _lambda_ufunc = lambda x: x
    else:
        print(f'This {norm_type} is not available for {name}')       
        
    
    return xr.apply_ufunc(_lambda_ufunc, da)#, input_core_dims=[['time']])

def regression(y, X, reg_names, units, config):
    nr = len(reg_names)
    ds_out = xr.apply_ufunc(_fit_ufunc, y, X, \
        input_core_dims = [['time'],['variable','time']], \
        output_core_dims=[['var','reg']], dask='allowed', \
        output_dtypes=[np.float], vectorize = True, \
        output_sizes = dict(var = 4, reg = nr), \
        kwargs = dict(rho = config['analysis']['rho']))

    #output processing
    ds_out = process_out(ds_out, reg_names[::-1], units, config)   
    return ds_out
    


# + {"code_folding": []}
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", default="config.yaml", help="Config yaml file")
    args = parser.parse_args()
    config_file = args.config
    print(config_file)
    config = yaml.load(open(config_file))
    if not os.path.exists(config["output_config"]["folder"]):
        os.mkdir(config["output_config"]["folder"])

    input_config = config['input_config']
    root_path = Path(input_config["folder"])
    s_year = config['analysis']['s_year']
    e_year = config['analysis']['e_year']
    sel_time_dict = dict(time = slice(str(s_year), str(e_year)))

    da_reg, reg_names = load_regressors(config, sel_time_dict)
    nr = len(reg_names)

    delayed_results = []
    for var in input_config['var']:
        for model in input_config['model']:
            for sim in input_config['simulation']:
                for ens in input_config['ens']:
                    # data opening
                    infile = list(root_path.glob(f'{var}_monthly_{model}_{sim}_{ens}_*.nc'))[0]
                    da_in = xr.open_dataset(infile)[var].squeeze(drop = True)#.chunk()
                    units = da_in.attrs['units']
                    da_in, da_in_time, da_in_month = process_in(da_in, var, sel_time_dict)

                    # regressors' postprocessing
                    da_reg['time'] = da_in_time
                    da_reg['month'] = da_in_month

                    # groupby splitting
                    analysis_type_ls = config['analysis']['type']
                    for analysis_type in analysis_type_ls:
                        if analysis_type == 'monthly':
                            y, X = da_in.groupby('time.month'), da_reg.groupby('time.month')
                        elif analysis_type == 'seasonally':
                            y, X = da_in.groupby('time.season'), da_reg.groupby('time.season')
                        elif analysis_type == 'annually':
                            y, X = da_in.copy(), da_reg.copy()
                        else:
                            print(f'None of this signal option is available')
                        # regression analysis
                        ds_out = dask.delayed(regression)(y, X, reg_names, units, config)
                        outfile = Path(config['output_config']['folder']) / f'{var}_{analysis_type}_{model}_{sim}_{ens}_{s_year}-{e_year}_results.nc'
                        delayed_obj = ds_out.to_netcdf(outfile, compute = False)
                        delayed_results.append(delayed_obj)

    with ProgressBar():
        results = dask.compute(*delayed_results)
    
if __name__ == "__main__":
    main()
