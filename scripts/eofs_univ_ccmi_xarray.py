import matplotlib
matplotlib.use('Agg')
import numpy as np
import time 
import supp_functions as fce
import matplotlib.pyplot as plt
from eofs.xarray import Eof
import statsmodels.api as sm
import sys
import re
import xarray as xr
import pandas as pd
import argparse
import platform
import os

pripona_nc = '.nc'

def xr_regression(y):
    X = sm.add_constant(reg, prepend=True) # regressor matrix
    mod = sm.GLSAR(y.values.squeeze(), X, 0, missing = 'drop') # MLR analysis without AR2 modeling
    res = mod.iterative_fit()

    return xr.DataArray(res.wresid)

def main(args):
    #environmental constants  
    if platform.system() == 'Windows':
        in_dir='../examples/'
        out_dir='../regressors/'
        reg_dir='../regressors/'#${in_dir}'regresory_2013/'
        nc_gen=True
        pdf_gen=False
        plus = ''
    else:
        n_samples = int(os.environ['n_samples'])
        in_dir = os.environ['in_dir']
        #out_dir = os.environ['out_dir']
        reg_dir = os.environ['reg_dir']
        pdf_gen = os.environ['pdf_gen']
        nc_gen = os.environ['nc_gen']

    what_re = args.what_re
    vari = args.vari
    i_year = args.i_year
    s_year = args.s_year
    e_year = args.e_year
    in_file_name = args.in_file_name

    if args.verbose:
        print('dataset: ', what_re)
        print('variable: ', vari)
        print('initial year of dataset: ', i_year)
        print('initial year of analysis: ', s_year)
        print('end year of analysis: ', e_year)
        print('input filename: ', in_file_name)


    print('data opening')
    in_netcdf = in_dir + in_file_name
    print(in_netcdf)
    ds = xr.open_dataset(in_netcdf)
    #print(ds)
    lat_name = fce.get_coords_name(ds, 'latitude')
    lat = ds.coords[lat_name]
    nlat = lat.shape[0]

    lev_name = fce.get_coords_name(ds, 'pressure')
    if ds.coords[lev_name].attrs['units'] == 'Pa':
        lev =  ds.coords[lev_name]/100.
        ds[lev_name] = lev    
    else:
        lev = ds.coords[lev_name]

    
    n = ds.coords['time'].shape[0]
    #it may happen that the field is 3D (longitude is missing)
    try:
        lon_name = fce.get_coords_name(ds, 'longitude')
        lon = ds.coords[lon_name]
        nlon = lon.shape[0]
    except:
        nlon = 1

    #print nlat, nlev, n, nlon

    #zonal mean
    if nlon != 1:
        uwnd = ds[vari].mean(lon_name)
    else:
        uwnd = ds[vari]
      
    #equatorial average and level selection
    sel_dict = {lev_name: fce.coord_Between(lev,10,50), lat_name: fce.coord_Between(lat,-10,10)}    
    zm_u = uwnd.sel(**sel_dict).mean(lat_name)
    #period selection
    times = pd.date_range(str(s_year)+'-01-01', str(e_year)+'-12-31', name='time', freq = 'M')
    zm_u_sel = zm_u.sel(time = times, method='ffill') #nearest
    #remove seasonality
    climatology = zm_u_sel.groupby('time.month').mean('time')
    anomalies = zm_u_sel.groupby('time.month') - climatology

    #print anomalies
    #sys.exit()
    
    #additional constants
    npca = 30
    norm=2 #5
    norms=3 #5
    what_sp = '' # what solar proxy?

    print("regressors' openning")
    global reg#, reg_names, nr
    reg, reg_names, history = fce.configuration_ccmi(what_re, what_sp, norm, 'no_qbo' , i_year, s_year, e_year, reg_dir)
    nr = reg.shape[1]
    #print(anomalies) 
    #extracting of other variability by MLR
    stacked = anomalies.stack(allpoints = [lev_name])
    stacked = stacked.reset_coords(drop=True)
    resids = stacked.groupby('allpoints').apply(xr_regression)
    resids = resids.rename({'dim_0': 'time'})
    resids['time'] = times
    #EOF analysis
    solver = Eof(resids.T, weights=None) 
    #sys.exit()

    #coslat = np.cos(np.deg2rad(lat)).clip(0.,1.)
    #wgts = np.sqrt(coslat)[np.newaxis,...]   

    for i in range(npca):		
            var_eofs = solver.varianceFraction(neigs=i)
            #print var_eofs
            if np.sum(var_eofs) > 0.95:
                    npca=i
                    total_variance = np.sum(var_eofs)
                    print(total_variance, ' % based on ', i, ' components')
                    break

    var_eofs = solver.varianceFraction(neigs=npca)
    pcs = solver.pcs(npcs=npca, pcscaling=1)
    nte = solver.northTest(neigs=npca, vfscaled=True)

    subdir = './'
    if pdf_gen:
        fig = plt.figure(figsize=(11,8))
        ax1 = fig.add_subplot(111)
        ax1.set_title(str(npca)+' PCAs cover '+str(np.round(total_variance*100, 2))+'% of total variance')
        for i in xrange(npca):
                #plotting
                pcs[:,i].plot(linewidth = 2, ax = ax1, label = 'pca '+str(i+1))

        ax1.set_xlabel('time [years]')
        ax1.set_ylabel('QBO index')
        ax1.set_title('')
        ax1.legend(loc = 'best')
        plt.savefig(reg_dir+'qbo_'+what_re+'_pcas.pdf', bbox_inches='tight')
        plt.close(fig)    
                
    if nc_gen:
        #save to netcdf
        #print(pcs[:,0])
        for i in range(npca):
            pcs_ds = pcs[:,i].to_dataset(name = 'index')
            pcs_ds.to_netcdf(reg_dir+r'qbo_'+what_re+'_pc'+str(i+1)+pripona_nc)       

if __name__ == "__main__":
    #inputs
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase output verbosity", action="store_true")
    parser.add_argument("what_re", help="dataset shor name")
    parser.add_argument("vari", help="variable name")
    parser.add_argument("i_year", help="initial year of dataset", type=int)
    parser.add_argument("s_year", help="initial year of analysis", type=int)
    parser.add_argument("e_year", help="end year of analysis", type=int)
    parser.add_argument("in_file_name", help="input filename")

    args = parser.parse_args()
    main(args)
    print('done')






