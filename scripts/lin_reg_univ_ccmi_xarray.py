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

def fit(y, X, reg_names):
    nr = len(reg_names)
    
    try:
        mod = sm.GLSAR(y.values, X, 2, missing = 'drop') # MLR analysis with AR2 modeling
        res = mod.iterative_fit()
        output = xr.Dataset({'coef': (['reg_name'], res.params[1:]), \
                'conf_int': (['reg_name', 'limit'], res.conf_int()[1:,:]), \
                'p_value': (['reg_name'],  res.pvalues[1:]), \
                'DWT': (sms.durbin_watson(res.wresid)), \
                'CoD': (res.rsquared)}, \
                coords = {'reg_name': (['reg_name'], reg_names),\
                          'limit': (['limit'], ['lower', 'upper'])})
    except: 
        nans = np.full([nr], np.nan)
        output = xr.Dataset({'coef': (['reg_name'], nans), \
                'conf_int': (['reg_name', 'limit'], np.array([nans, nans]).T), \
                'p_value': (['reg_name'],  nans), \
                'DWT': (np.nan), \
                'CoD': (np.nan)}, \
                coords = {'reg_name': (['reg_name'], reg_names),\
                          'limit': (['limit'], ['lower', 'upper'])})

    return output

def xr_regression(y, **kwargs):
    X = sm.add_constant(np.array(kwargs['reg']), prepend=True) # regressor matrix
    res_ls = []
    nr = kwargs['nr']
    n = y.shape[0]
    datum = range(1,13)
    datum *= (n/12)
    if kwargs['monthly']:
        n_iter = 13
    else:
        n_iter = 1

    for mi in xrange(n_iter):
        monthi = mi == np.array(datum)
        if mi == 0:
            res = fit(y[~monthi], X[~monthi], kwargs['reg_names'])
        else:
            res = fit(y[monthi], X[monthi], kwargs['reg_names'])
        res_ls.append(res)
    res_da = xr.concat(res_ls, dim = 'month')
    res_da['month'] = np.arange(0, n_iter)
    return res_da

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
    monthly = args.monthly
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
        conf_str = conf_str[:-3] # reg. conf. is the same for filt. or unfilt. analysis
        filt_years = [1982,1983,1984]
    elif conf_str[-2:] == 'pi':
        conf_str = conf_str[:-3]
        filt_years = [1991,1992,1993]
    elif conf_str[-2:] == 'bo':
        conf_str = conf_str[:-3]
        filt_years = [1982,1983,1984,1991,1992,1993]
    else:
        filt_years = None
        
    print('data opening')
    in_netcdf = in_dir + in_file_name 
    #ds = xr.open_mfdataset(in_netcdf, concat_dim = 'ens')
    ds = xr.open_dataset(in_netcdf)
    #ds = ds.chunk({'lat': 10000})
    lat_name = fce.get_coords_name(ds, 'latitude')
    lat = ds.coords[lat_name].values
    nlat = lat.shape[0]


    lev_name = fce.get_coords_name(ds, 'air_pressure')
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

    #currently is tested for zonally averaged files only
    if nlon != 1:
        ds = ds.mean(lon_name)

    print("regressors' openning")
    #global reg, reg_names, nr
    reg, reg_names, history = fce.configuration_ccmi(what_re, what_sp, norm, conf_str, i_year, s_year, e_year, reg_dir, filt_years = filt_years)
    nr = reg.shape[1]
    #select date range and variable
    #times = pd.date_range(str(s_year)+'-01-01', str(e_year)+'-12-31', name='time', freq = 'M')
    ds_sel = fce.date_range_xr(ds, i_year, s_year, e_year, 1, 12, n, filt_years = filt_years)#.sel(lat = slice(-25,25))#ds.sel(time = times, method='ffill') #nearest #[vari]
    print('anomalies calculation')
    anomalies, _ = fce.deseasonalize(ds_sel)
    anomalies = anomalies.squeeze().reset_coords(drop=True)
    print('regression calculation')
    filt_dims =  list(filter(lambda x: x not in ['time'], anomalies.dims))[::-1]
    stacked = anomalies[vari].stack(allpoints = filt_dims)#[lev_name, lat_name])

    reg_kwargs = dict(reg = reg, reg_names = reg_names, nr = nr, monthly = monthly)
    coefs = stacked.groupby('allpoints').apply(xr_regression, **reg_kwargs)#.squeeze()

    if lev_name in anomalies.dims:
        coefs['allpoints']  = stacked.coords['allpoints'].sortby(lev_name) # I need to sort allpoints multiindex according to lev, otherwise I would get reversedcoefs   
    else:
        coefs['allpoints']  = stacked.coords['allpoints']
    
    ndims = len(filt_dims)
    if ndims != 1:
        cu_ds = coefs.unstack('allpoints')
    else:
        cu_ds = coefs.rename({'allpoints': filt_dims[0]})
    
    if fce.str2bool(nc_gen):
        print('netCDF Output')
        cu_ds['coef'].attrs['long_name'] = 'Regression coefficients'
        cu_ds['CoD'].attrs['long_name'] = 'Coefficient of determination'
        cu_ds['p_value'].attrs['long_name'] = 'Statistical significance (p-value)'
        cu_ds['DWT'].attrs['long_name'] = 'Durbin-Watson test'
        if  filt_dims in [lev_name]:
            cu_ds[lev_name].attrs['units'] = 'hPa'
        cu_ds.attrs['history'] = 'Regressors included: '+history
        cu_ds.attrs['description'] = 'Created ' + time.ctime(time.time()) 
        cu_ds.to_netcdf(out_dir+'stat_outputs'+suffix_nc)
    
    if fce.str2bool(pdf_gen) and set(filt_dims) == set([lat_name, lev_name]):
        print('RC visualization')
        my_cmap = mpl.colors.ListedColormap(['yellow', 'red', 'white'])
        fgp = xr.plot.FacetGrid(cu_ds['p_value'], row = 'reg_name', col = 'month', sharey = True, sharex = True)
        plot_cf_kwargs = dict(yincrease=False, levels = [0,0.01,0.05,1], add_colorbar = False, cmap=my_cmap)
        fgp.map_dataarray(xr.plot.contourf, lat_name, lev_name, **plot_cf_kwargs) 
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

        plot_kwargs_zero = dict(yincrease=False, cmap=('k'), linewidths = 6, add_colorbar=False, levels = [0])
        plot_kwargs = dict(yincrease=False, cmap=('k'), add_colorbar=False, levels = c_levels[c_levels>0], linewidths = 3)
        fgp.data = cu_ds['coef']
        fgp.map_dataarray(xr.plot.contour, lat_name, lev_name, **plot_kwargs_zero)
        
        fgp.data = cu_ds['coef']
        fgp.map_dataarray(xr.plot.contour, lat_name, lev_name, **plot_kwargs)

        plot_kwargs['levels'] = c_levels[c_levels<0]
        plot_kwargs['linestyles'] = 'dashed'
        fgp.data = cu_ds['coef']
        fgp.map_dataarray(xr.plot.contour, lat_name, lev_name, **plot_kwargs)

        ax = fgp.axes[0,0]                                                          
        #ax.set_ylabel('pressure [hPa')
        #ax.set_xlabel('latitude [deg]')
        ax.set_ylim(1000,0.1)
        ax.set_yscale('log')

        plt.savefig(out_dir+'visualization'+suffix_pdf, bbox_inches = 'tight')

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
    choices_def = ['all_trend','all_2trends','all_eesc', 'no_saod_trend', 'no_saod_2trends', 'no_saod_eesc', 'no_saod_enso_trend', 'no_saod_enso_2trends', 'no_saod_enso_eesc', 'massi_trend', 'massi_2trends', 'massi_notrend', 'massi_eesc', 'solarAp_trend']
    choices = choices_def + [i+'_bo' for i in choices_def] + [i+'_pi' for i in choices_def] + [i+'_el' for i in choices_def]
    parser.add_argument("config", help="regression configuration", choices=choices)
    parser.add_argument("--monthly", dest = 'monthly', action = 'store_true')
    #parser.set_defaults(monthly = False)
    #parser.add_argument("zonal_config", help="Zonal, whole or map?")   


    args = parser.parse_args()
    main(args)
    print('{} seconds elapsed'.format(time.time()-start))
    print('done')



