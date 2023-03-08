import xarray as xr
import panel as pn
import panel.widgets as pnw
import numpy as np

def watl_anim(da, quant=0.9):
    """ Plot time series of water level contours.
    
    Args:
        da(xr DataArray): DataArray with water levels (in m).
        quant (float, optional): maximum absolute colour bar level (in % of absolute water levels). Default to 0.9.
    Returns:
        panel object with water level animation.
    """
    vmax=(xr.apply_ufunc(np.abs,da)).quantile(quant)
    t = pnw.Player(name='t',start=0,end=len(da.t)-1,step=1,loop_policy='loop')
    fig = lambda t:da.isel(t=t).hvplot(cmap="Blues",clim=(-vmax,vmax))
    a=pn.interact(fig,t=t)
    return pn.Column(pn.Row(a[1]),pn.Row(a[0]))