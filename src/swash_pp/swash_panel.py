import xarray as xr
import panel as pn
import panel.widgets as pnw
import numpy as np
import math

def watl_anim(ds,quant=0.9,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=1,vel_clip_max=0.9,vel_clip_min=0.1,):
    """ Plot time series of water level contours.
    
    Args:
        ds (xr Dataset): Dataset with 'Watlev' (in m), 'Vksi', and 'Veta' (in m/s).
        quant (float, optional): maximum absolute colour bar level (in % quantile of absolute water levels). Default to 0.9.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        vel_clip_max (float, optional): maximum x- and y-velocity clip (in % of mag quantile). Default to 0.9.
        vel_clip_min (float, optional): minimum absolute velocity (in % of mag quantile) to be plotted (otherwise nan). Default to 0.1          
    
    Returns:
        panel object with water level animation.
    """
    # Assign xmin, xmax, dx, ymin, ymax, and dy if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    if "Vksi" in ds:
        # define mag and angle - input for bokeh vectorfield
        temp["mag"]=(temp["Vksi"]**2+temp["Veta"]**2)**0.5
        temp["angle"]=xr.apply_ufunc(np.arctan2,temp["Veta"],temp["Vksi"])
        # velocity clipping
        vel_clip_max=xr.apply_ufunc(np.abs,temp["mag"]).quantile(vel_clip_max).item()
        temp["mag"]=temp["mag"].clip(max=vel_clip_max)
        temp=temp.where(lambda x:x["mag"]>=vel_clip_min,drop=True)
        temp["mag"]=temp["mag"]*scale
    # water level colour bar
    vmax=(xr.apply_ufunc(np.abs,ds.Watlev)).quantile(quant) # absolute max colourbar level
    # panel widget
    t = pnw.Player(name='t',start=0,end=len(ds.t)-1,step=1,loop_policy='loop')
    if "Vksi" in ds:
        fig = lambda t:ds.Watlev.sel(x=slice(xmin,xmax),y=slice(ymin,ymax)).isel(t=t).hvplot(cmap="Blues",clim=(-vmax,vmax),clabel='Water level [m]').opts(colorbar_opts={'title_standoff': -150,'padding': 30}) *\
                        temp.isel(t=t).hvplot.vectorfield(x='x', y='y', angle='angle', mag='mag', hover=False).opts(magnitude='mag')
    else:
        fig = lambda t:ds.Watlev.sel(x=slice(xmin,xmax),y=slice(ymin,ymax)).isel(t=t).hvplot(cmap="Blues",clim=(-vmax,vmax),clabel='Water level [m]').opts(colorbar_opts={'title_standoff': -150,'padding': 30})
    a=pn.interact(fig,t=t)
    return pn.Column(pn.Row(a[1]),pn.Row(a[0]))