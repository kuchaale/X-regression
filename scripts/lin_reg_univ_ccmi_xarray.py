import matplotlib                                                       
matplotlib.use('Agg')
import numpy as np
import supp_functions as fce
import sys
import os
import xarray as xr
import pandas as pd
import argparse
import time
import platform
import statsmodels.api as sm
import statsmodels.stats.stattools as sms
import matplotlib as mpl
import matplotlib.pyplot as plt


periods=['','_01_jan','_02_feb','_03_mar','_04_apr','_05_may','_06_jun','_07_jul','_08_aug','_09_sep','_10_oct','_11_nov','_12_dec','_win', '_spr', '_sum', '_aut']
bool_str = ['true', '1', 't', 'y', 'yes', 'yeah', 'yup', 'certainly', 'uh-huh']

def xr_regression(y):
    X = sm.add_constant(reg, prepend=True) # regressor matrix
    #nr = reg.shape[1]
    #print(y)
    try:
        mod = sm.GLSAR(y.values, X, 2, missing = 'drop') # MLR analysis with AR2 modeling
        res = mod.iterative_fit() 
        output = xr.DataArray([res.params[1:], res.pvalues[1:], [sms.durbin_watson(res.wresid)]*nr, [res.rsquared]*nr], coords=[['coefs', 'p_values', 'dwt', 'cod'], reg_names], dims=['stat_var', 'regs'])
    except: 
        nans = np.full([nr], np.nan)
        output = xr.DataArray([nans, nans, nans, nans], coords=[['coefs', 'p_values', 'dwt', 'cod'], reg_names], dims=['stat_var', 'regs'])

    return output

def main(args):
    #environmental constants  
    if platform.system() == 'Windows':
        n_samples=5000
        in_dir='../examples/'
        out_dir=''
        reg_dir='../regressors/'#${in_dir}'regresory_2013/'
        pdf_gen=True
        nc_gen=True
    else:
        n_samples = int(os.environ['n_samples'])
        in_dir = os.environ['in_dir']
        out_dir = os.environ['out_dir']
        reg_dir = os.environ['reg_dir']
        pdf_gen = os.environ['pdf_gen']
        nc_gen = os.environ['nc_gen']

    plus = ''
    # additional constants
    what_sp = '' # what solar proxy?
    norm = 4 # normalization type
    

    """Run the program."""
    what_re = args.what_re
    vari = args.vari
    i_year = args.i_year
    s_year = args.s_year
    e_year = args.e_year
    in_file_name = args.in_file_name
    conf_str = args.config
    suffix_pdf = '_{}-{}_{}.pdf'.format(s_year, e_year, conf_str)
    suffix_nc = '_{}-{}_{}.nc'.format(s_year, e_year, conf_str)
    
    #zonal_b = args.zonal_config
    

    out_dir += vari+'_'+what_re+'_'

    if args.verbose:
        print('dataset: ', what_re)
        print('variable: ', vari)
        print('initial year of dataset: ', i_year)
        print('initial year of analysis: ', s_year)
        print('end year of analysis: ', e_year)
        print('input filename: ', in_file_name)
        print('regression configuration: ', conf_str)
        #print('Zonal, whole or map? ', zonal_b)
   

    if conf_str[-2:] == 'el':
        filt_years = [1982,1983,1984]
    elif conf_str[-2:] == 'pi':
        filt_years = [1991,1992,1993]
    elif conf_str[-2:] == 'bo':
        filt_years = [1982,1983,1984,1991,1992,1993]
    else:
        filt_years = None
        
    print('data opening')
    in_netcdf = in_dir + in_file_name 
    #ds = xr.open_mfdataset(in_netcdf, concat_dim = 'ens')
    ds = xr.open_dataset(in_netcdf)
    
    lat_name = fce.get_coords_name(ds, 'latitude')
    lat = ds.coords[lat_name].values
    nlat = lat.shape[0]


    lev_name = fce.get_coords_name(ds, 'pressure')
    if ds.coords[lev_name].attrs['units'] == 'Pa':
        lev =  ds.coords[lev_name].values/100.
        ds[lev_name] = lev    
    else:
        lev = ds.coords[lev_name].values
    #print(lev) 
    nlev = lev.shape[0]
    gen = np.arange(nlev)

    n = ds.coords['time'].shape[0]

    #it may happen that the field is 3D (longitude is missing)
    try:
        lon_name = fce.get_coords_name(ds, 'longitude')
        lon = ds.coords[lon_name].values
        nlon = lon.shape[0]
    except:
        nlon = 1

    if nlon != 1:
        ds = ds.mean(lon_name)

    print("regressors' openning")
    global reg, reg_names, nr
    reg, reg_names, history = fce.configuration_ccmi(what_re, what_sp, norm, conf_str[:-3], i_year, s_year, e_year, filt_years = filt_years)
    X = reg[:,:] 
    nr = X.shape[1]

    #select date range and variable
    #times = pd.date_range(str(s_year)+'-01-01', str(e_year)+'-12-31', name='time', freq = 'M')
    ds_sel = fce.date_range_xr(ds, i_year, s_year, e_year, 1, 12, n, filt_years = filt_years)#ds.sel(time = times, method='ffill') #nearest #[vari]

    print('anomalies calculation')
    anomalies, _ = fce.deseasonalize(ds_sel)
    anomalies = anomalies.squeeze().reset_coords(drop=True)
    print('regression calculation')
    #if 'ens' in ds.dims.keys():
    #    stacked = anomalies.stack(allpoints = ['lev', 'lat', 'lon', 'ens'])
    #else:
    #    stacked = anomalies.stack(allpoints = ['lev', 'lat', 'lon'])
    stacked = anomalies[vari].stack(allpoints = [lev_name, lat_name])
    #print(anomalies)
    #stacked = anomalies[vari].stack(allpoints = anomalies.dims.keys()[:-1])

    #stacked = stacked.reset_coords(drop=True)
    coefs = stacked.groupby('allpoints').apply(xr_regression)
    coefs['allpoints']  = stacked.coords['allpoints'].sortby(lev_name) # I need to sort allpoints multiindex according to lev, otherwise I would get reversedcoefs   
    coefs_unstacked = coefs.unstack('allpoints')
    print('output processing')
    cu_ds = coefs_unstacked.to_dataset(dim = 'stat_var')

    cod = cu_ds.cod.isel(regs=[0]).squeeze()
    dwt = cu_ds.dwt.isel(regs=[0]).squeeze()

    cu_ds = fce.subset_variables(cu_ds, vlist = ['p_values', 'coefs'])

    cod_ds = cod.to_dataset(name = 'cod')
    cod_ds.reset_coords(drop=True, inplace=True)
    dwt_ds = cod.to_dataset(name = 'dwt')
    dwt_ds.reset_coords(drop=True, inplace=True)

    cu_ds = cu_ds.merge(cod_ds)
    cu_ds = cu_ds.merge(dwt_ds)
    
    if nc_gen:
        print('netCDF Output')
        cu_ds.coefs.attrs['long_name'] = 'Regression coefficients'
        cu_ds.cod.attrs['long_name'] = 'Coefficient of determination'
        cu_ds.p_values.attrs['long_name'] = 'Statistical significance (p-value)'
        cu_ds.dwt.attrs['long_name'] = 'Durbin-Watson test'
        #cu_ds.coords[lev_name].attr['units'] = 'hPa'
        cu_ds.attrs['history'] = 'Regressors included: '+history
        cu_ds.attrs['description'] = 'Created ' + time.ctime(time.time()) 
        cu_ds.to_netcdf(out_dir+'stat_outputs'+suffix_nc)
    

    if pdf_gen:
        print('solar RC visualization')
        fig, ax = plt.subplots(figsize=(12,9))
        my_cmap = mpl.colors.ListedColormap(['yellow', 'red', 'white'])
        coefs_unstacked.sel(stat_var = 'p_values', regs = 'solar').squeeze().plot.contourf(yincrease=False, levels = [0,0.01,0.05], cmap=my_cmap, ax = ax)
        if vari in ['zmta']:
            c_levels = [-30,-15,-10,-5,-2,-1,-0.5,-0.25]
            c_levels += [0]+fce.rev_sign(c_levels)
            c_levels = np.array(c_levels)
        elif vari in ['zmua']:
            c_levels = [-30,-15,-10,-5,-2,-1]
            c_levels += [0]+fce.rev_sign(c_levels)
            c_levels = np.array(c_levels)
        else:
            c_levels = np.arange(-10,11,1)

        plot_kwargs_zero = dict(yincrease=False, cmap=('k'), linewidths = 6, add_colorbar=False, levels = [0], ax = ax)
        plot_kwargs = dict(yincrease=False, colors='k', add_colorbar=False, levels = c_levels[c_levels>0], ax = ax, linewidths = 3)
        coefs_unstacked.sel(stat_var = 'coefs', regs = 'solar').squeeze().plot.contour(**plot_kwargs_zero)
        coefs_unstacked.sel(stat_var = 'coefs', regs = 'solar').squeeze().plot.contour(**plot_kwargs)
        plot_kwargs['levels'] = c_levels[c_levels<0]
        plot_kwargs['linestyles'] = 'dashed'
        coefs_unstacked.sel(stat_var = 'coefs', regs = 'solar').squeeze().plot.contour(**plot_kwargs)
        ax.set_yscale('log')
        ax.set_ylabel('pressure [hPa')
        ax.set_xlabel('latitude [deg]')
        ax.set_ylim(1000,0.1)
        plt.savefig(out_dir+'visualization'+suffix_pdf, bbox_inches = 'tight')
        plt.close(fig)

if __name__ == "__main__":
    start = time.time()
    #inputs
    description = 'MLR input arguments.'
    parser = argparse.ArgumentParser(description = description)
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("what_re", help="dataset shor name")
    parser.add_argument("vari", help="variable name")
    parser.add_argument("i_year", help="initial year of dataset", type=int)
    parser.add_argument("s_year", help="initial year of analysis", type=int)
    parser.add_argument("e_year", help="end year of analysi", type=int)
    parser.add_argument("in_file_name", help="input filename")
    choices_def = ['all_trend','all_2trends','all_eesc', 'no_saod_trend', 'no_saod_2trends', 'no_saod_eesc', 'no_saod_enso_trend', 'no_saod_enso_2trends', 'no_saod_enso_eesc']
    choices = choices_def + [i+'_bo' for i in choices_def] + [i+'_pi' for i in choices_def] + [i+'_el' for i in choices_def]
    parser.add_argument("config", help="regression configuration", choices=choices)
    #parser.add_argument("zonal_config", help="Zonal, whole or map?")   


    args = parser.parse_args()
    main(args)
    print('{} seconds elapsed'.format(time.time()-start))
    print('done')



