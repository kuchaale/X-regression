from netCDF4 import Dataset
import numpy as np
from time import time, ctime
import os
import calendar
import statsmodels.api as sm


def deseasonalize(ds, how = 'absolute'):
    time_var = get_coords_name(ds, 'time')
    climatology = ds.groupby(time_var+'.month').mean(time_var)#.median(time_var) does not work
    anomalies = ds.copy()
    for key in ds.data_vars.keys():
        if ds[key].ndim > 2:
            if key.lower() in ['vmro3', 'o3']:
                anomalies[key] = (ds[key].groupby(time_var+'.month') - climatology[key]).groupby(time_var+'.month')/climatology[key]*100
            else:
                anomalies[key] = (ds[key].groupby(time_var+'.month') - climatology[key])

    return anomalies, climatology

def subset_variables(ds, vlist = ['TMP', 'VGRD', 'UGRD']):     
    """
    Reduces an xarray dataset ds to only contain the variables in vlist.

    Parameters
    ----------
    ds : xarray Dataset

    Returns
    -------
    ds : xarray Dataset with subset variables

    Notes
    -----

    Author
    ------
    Phillip J. Wolfram
    05/05/2016

    Examples
    --------
    """

    # get set of variables to drop (all ds variables not in vlist)
    dropvars = set(ds.data_vars.keys()) - set(vlist)
    
    # drop spurious variables
    ds = ds.drop(dropvars)

    # must also drop all coordinates that are not associated with the variables
    coords = set()
    for avar in ds.data_vars.keys():
        coords |= set(ds[avar].coords.keys())
    dropcoords = set(ds.coords.keys()) - coords

    # drop spurious coordinates
    ds = ds.drop(dropcoords)
    
    return ds


def linreg(y, X):
	mod = sm.GLSAR(y, X, 2, missing = 'drop')
	res = mod.iterative_fit()
	return res.params[1:], res.pvalues[1:], np.array(res.conf_int())[1:,:], res.rsquared, res.fittedvalues, res.wresid

def my_function(y, X):
    return linreg(y, X)[0]

def get_var(ds, l_name):
    c_name = [k for k, v in ds.data_vars.iteritems() if 'long_name' in v.attrs.keys() and l_name in v.long_name][0]
    return ds.data_vars[c_name].values
def get_var_name(ds, l_name):
    return [k for k, v in ds.data_vars.iteritems() if 'long_name' in v.attrs.keys() and l_name in v.long_name][0]

def get_coords(ds, l_name):
    c_name = [k for k, v in ds.coords.iteritems() if 'long_name' in v.attrs.keys() and l_name in v.long_name][0]
    #if l_name == 'pressure' and ds.coords[c_name].attrs['units'] == 'Pa':
    #    return ds.coords[c_name].values/100.
    #else:
    return ds.coords[c_name].values

def get_coords_name(ds, l_name):
    return [k for k, v in ds.coords.iteritems() if 'standard_name' in v.attrs.keys() and l_name in v.standard_name][0]

def configuration_ccmi(what_re, what_sp, norm, conf, i_year, s_year, e_year, filt_years = None):
    reg_dir = ''#os.environ['reg_dir']
    print reg_dir
    qbo1 = open_reg_ccmi(reg_dir+'qbo_'+what_re+what_sp+'_pc1.nc', 'index', norm, i_year, s_year, e_year, filt_years = filt_years)
    n = qbo1.shape[0]
    qbo2 = open_reg_ccmi(reg_dir+'qbo_'+what_re+what_sp+'_pc2.nc', 'index', norm, i_year, s_year, e_year, filt_years = filt_years)
    solar = open_reg_ccmi(reg_dir+'solar_1947.nc', 'solar', 0, 1947, s_year, e_year, filt_years = filt_years)
    solar /= 126.6 # normalization on Smax - Smin sfu
    what_re2 = 'HadISST'
    i_year2 = 1947
    e_year2 = 2009
    enso = open_reg_ccmi(reg_dir+'enso_'+what_re2+'_monthly_'+str(i_year2)+'_'+str(e_year2)+'.nc', 'enso', norm, i_year2, s_year, e_year, filt_years = filt_years)
    if what_re in ['20CR','era-40','ncep1','m_iau','era_interim_mplev']: #'jra55',
        saod = open_reg_ccmi(reg_dir+'saod_1960_2013.nc', 'saod', norm, 1960, s_year, e_year, filt_years = filt_years)
    else:
        saod = open_reg_ccmi(reg_dir+'sad_gm_50hPa_1949_2013.nc', 'sad', 0, 1949, s_year, e_year, filt_years = filt_years)
    #if vari in ['O3'] and s_year >= 1979 and conf_str[:7] != '2trends':
    #eesc = open_reg_ccmi(reg_dir+'eesc.nc', 'eesc', norm, 1979, s_year, e_year, filt_years = filt_years)
    trend = np.linspace(-1, 1, n)

    #elif vari == 'x':
    infl_point = 1996
    year=np.concatenate([[i_year+i]*12 for i in range((n)/12)])
    part1 =  (year <= infl_point)
    part2 =  (year >= infl_point)
        
    trend1 = np.zeros(n)
    trend2 = np.zeros(n)

    trend1[part2] = np.linspace(0,1,np.count_nonzero(part2))
    trend2[part1] = np.linspace(-1,0,np.count_nonzero(part1))

    volc_i = (year >= s_year) & (year <= e_year)
    trend1 =  trend1[volc_i]
    trend2 =  trend2[volc_i]     

    if conf == 'all_trend':                     
        reg = np.column_stack((trend, solar, enso, saod, qbo1, qbo2)) 
	my_xticks = ['trend', 'solar', 'enso', 'saod', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'all_2trends':
        reg = np.column_stack((trend1, solar, enso, saod, qbo1, qbo2, trend2)) 
	my_xticks = ['trend', 'solar', 'enso', 'saod', 'qbo1', 'qbo2', 'trend2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'all_eesc':
        reg = np.column_stack((eesc, solar, enso, saod, qbo1, qbo2)) 
	my_xticks = ['eesc', 'solar', 'enso', 'saod', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_trend':
        reg = np.column_stack((trend, solar, enso, qbo1, qbo2)) 
	my_xticks = ['trend', 'solar', 'enso', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_2trends':
        reg = np.column_stack((trend1, solar, enso, qbo1, qbo2, trend2)) 
	my_xticks = ['trend', 'solar', 'enso', 'qbo1', 'qbo2', 'trend2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_eesc':
        reg = np.column_stack((eesc, solar, enso, qbo1, qbo2)) 
	my_xticks = ['eesc', 'solar', 'enso', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_enso_trend':
        reg = np.column_stack((trend, solar, qbo1, qbo2)) 
	my_xticks = ['trend', 'solar', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_enso_2trends':
        reg = np.column_stack((trend1, solar, qbo1, qbo2, trend2)) 
	my_xticks = ['trend', 'solar', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)
    elif conf == 'no_saod_enso_eesc':
        reg = np.column_stack((eesc, solar, qbo1, qbo2)) 
	my_xticks = ['eesc', 'solar', 'qbo1', 'qbo2']
	#outdir += '/rel_impact/LIN_REG/'
	history =  ", ".join(my_xticks)

    return reg, my_xticks, history


def date_range_ccmi(arr, s_year_data, s_year, e_year, s_mon, e_mon):
    n = arr.shape[0]
    year=np.concatenate([[s_year_data+i]*12 for i in range(np.int(np.ceil(n/12.)))])
    month = range(1,13)*np.int(np.ceil(n/12.))
    #print year
    year = year[:n]
    month = month[:n]
    #print month
    #print year
    volc_i = (year >= s_year) & (year <= e_year)
    n_dim = np.array(arr).ndim
    #print type(s_year), type(e_year)
    #print year
    #print volc_i
    #print year[volc_i].shape
    #print type(arr), arr.shape
    if n_dim > 1:
        arr = np.array(arr[volc_i,...])
    else:
        arr = np.array(arr[volc_i])

    #print arr.shape
    index1 = s_mon-1
    if e_mon != 12 and month[-1] == 12:
        index2 = -(12-e_mon)
    else:
        index2 = year.shape[0]

    if n_dim > 1:
        res = arr[index1:index2,...]
    else:
        res = arr[index1:index2]

    return res

def date_range_xr(arr, s_year_data, n, s_year, e_year, s_mon, e_mon, filt_years = None):
    times = pd.date_range(start=str(s_year_data), periods=n, freq=pd.tseries.offsets.DateOffset(months=1))
    #print(arr)
    arr['time'] = times
    tr = pd.date_range(start=str(s_year)+'-'+str(s_mon), end = str(e_year)+'-'+str(e_mon), freq=pd.tseries.offsets.DateOffset(months=1))
    if filt_years != None:
        tr = filter(lambda x: x.year not in filt_years, tr)
    arr = arr.sel(time = pd.to_datetime(tr))
    return arr

def normalize(data, norm, down_bound = -1., upper_bound = 1.):
    avg = np.nanmean(data)
    sd = np.nanstd(data)
    if norm == 0:
        data = data
    elif norm == 1:
        data = (data-avg)/sd
    elif norm == 2:
        data = (data-avg)/avg
    elif norm == 3:
	#print avg
        data = data/avg
    elif norm == 4:
        data = data/sd
    elif norm == 5:
        #print data
        dh = np.nanmax(data)
        dl = np.nanmin(data)
        #print dh
        #print dl
        data = (((data-dl)*(upper_bound-down_bound))/(dh-dl))+down_bound
        #print data
    elif norm == 6:    
	data = data-avg
    elif norm == 7:	
	data = data - runningMeanFast(data,6)
    return np.array(data)

def normalize_xr(data, norm, down_bound = -1., upper_bound = 1.):
    avg = data.mean('time')
    sd = data.std('time')
    if norm == 0:
        data = data
    elif norm == 1:
        data = (data-avg)/sd                                                               
    elif norm == 2:                                                                        
        data = (data-avg)/avg
    elif norm == 3:
        data = data/avg
    elif norm == 4:
        data = data/sd
    elif norm == 5:                                                                        
        dh = data.max()
        dl = data.min()
        #print dl                                                                          
        data = (((data-dl)*(upper_bound-down_bound))/(dh-dl))+down_bound
        #print data
    elif norm == 6:                                                                        
        data = data-avg

    return data	

def open_reg_ccmi(in_file, var, norm, i_year, s_year, e_year, filt_years = None):
    ds = xr.open_dataset(in_file, decode_times=False)
    data_xr = normalize_xr(ds[var], norm)
    s_mon = 1
    e_mon = 12    
    data = date_range_xr(data, i_year, s_year, e_year, s_mon, e_mon, filt_years = filt_years)
    return data

if __name__ == '__main__':
    main()


