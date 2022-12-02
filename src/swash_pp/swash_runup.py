"""
Created on Nov 30 2022
Scripts to analyse runup.
@author: rfonsecadasilva
"""

def swash_R2(Ibp):
    """Calculate 2% runup from instantaneous water level using zero-downcrossing method.

    Args:
        Ibp (xr data array): Instantaneous (vertical) beach position (in m relative to stil water level).

    Returns:
        np array: 2% runup
    """
    fir=Ibp.isel(t=slice(None,-1))-Ibp.mean("t") # first
    sec=Ibp.isel(t=slice(1,None)).assign_coords(t=fir.t)-Ibp.mean("t")# second
    time_cross=fir.where((fir*sec<0) & (fir>0),drop=True) # calculate crossing points for zero-downcrossing on swash
    R2=Ibp.groupby_bins('t',time_cross.t).max().quantile(0.98).values #calculate 2% runup
    return R2