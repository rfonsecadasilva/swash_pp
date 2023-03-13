"""
Created on Nov 30 2022
Scripts to analyse runup.
@author: rfonsecadasilva
"""
import xarray as xr
import numpy as np

def swash_R2(Ibp,exc_value=None):
    """Calculate 2% runup from instantaneous water level using zero-downcrossing method.

    Args:
        Ibp (xr data array): Instantaneous (vertical) beach position (in m relative to stil water level).
        exc_value (float): Exception value for runup (abs values greater than exc_value are excluded; in m). Default to None.

    Returns:
        R2 (xr data array): xr data array with R2.
    """
    if exc_value:
        bad_data=(xr.apply_ufunc(np.isnan,Ibp.where(lambda x:np.abs(x)>exc_value)).sum()/Ibp.size*100).values.item()
        if bad_data!=100:
            print(f'{100-bad_data:.1f} % of runup data is corrupted')
        Ibp=Ibp.where(lambda x:np.abs(x)<exc_value,drop=True)
    time_cross=Ibp.where(( (Ibp.shift(t=1) - Ibp.mean(dim="t")) * (Ibp - Ibp.mean(dim="t")) <0) & (Ibp - Ibp.mean(dim="t") >0),drop=True) # calculate crossing points for zero-downcrossing on swash
    if Ibp.y.size==1:
        if "y" in Ibp.dims: Ibp=Ibp.isel(y=0)
        R2=Ibp.groupby_bins('t',time_cross.t).max(dim="t").quantile(0.98,dim="t_bins").values #calculate 2% runup
        R2=xr.DataArray(R2[0], dims=())
    else:
        R2=[]
        for i in range(len(Ibp.y)): #iteration is needed becase groupby_bins method cannot consider 2D groups to perform this operation
            R2.append(Ibp.isel(y=i).groupby_bins('t',time_cross.isel(y=i).dropna(dim="t").t).max(dim="t").quantile(0.98,dim="t_bins").values.item())
        R2=xr.DataArray(R2, dims=("y"))
        R2=R2.assign_coords(y=Ibp.y)
    R2.attrs["standard_name"],R2.attrs["long_name"],R2.attrs["units"]="2% runup","2% runup","m"
    return R2