import xarray as xr
from scipy.signal import butter, sosfilt

def filter_lf_hf(da,fcut):
    """ Filter time series into high and low frequency bands.
    Note that time averaged of time series is removed.

    Args:
        da (xr DataArray): DataArray with variable to be filtered.
        fcut (float): cutoff frequency (in Hz)

    Returns:
        da_lf (xr Datarray): DataArray with low-frequency variable.
        da_hf (xr Datarray): DataArray with high-frequency variable.
    """
    fs=1/(da.t.isel(t=1)-da.t.isel(t=0)).item() # output frequency
    # define low- and high-pass 5th order Butterworth filter
    sos_lf=butter(5, fcut, 'lp', fs=fs,output='sos')
    sos_hf=butter(5, fcut, 'hp', fs=fs,output='sos')    
    def butter_filter(data, sos): # create a ufunc to apply the filter
        return sosfilt(sos, data)
    da_lf=xr.apply_ufunc(butter_filter,da-da.mean("t"),sos_lf)
    da_hf=xr.apply_ufunc(butter_filter,da-da.mean("t"),sos_hf)
    return da_lf,da_hf
