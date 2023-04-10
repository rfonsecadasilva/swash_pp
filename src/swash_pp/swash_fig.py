import warnings, math
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import swash_pp.lwt as lwt
import string

def Hs_Mfv_fig(ds,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=None,vel_clip_max=0.9,vel_clip_min=0.1,plot_dep_dev=False,dep_levels=None,hs_levels=None,hs_ticks=None,Hs0=None,cmap="jet"):
    """
    Create 2D fig with significant wave height and mass flux velocities.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'Hs', 'Mfvx' and Mfvy'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        vel_clip_max (float, optional): maximum x- and y-velocity clip (in % of x quantile). Default to 0.9.
        vel_clip_min (float, optional): minimum absolute velocity (in % of x quantile) to be plotted (otherwise nan). Default to 0.1
        plot_dep_dev (bool), optional: Condition plotting deviations from depth at the first cross-shore section. Default to False.        
        dep_levels(np array, optional): array with depth contour levels (in m) to be plotted. If None, no depth contours.
        hs_levels(np array, optional): array with significant wave height contour levels (in m) to be plotted.
        hs_ticks(np array, optional): array with significant wave height contour levels ticks (in m) to be plotted.
        Hs0 (float, optional): deep water significant wave height (in m). If provided, Hs is normalised.
        cmap (str): hs matplotlib colour map. Default to "jet".
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, dx, ymin, ymax, dy, tmin, tmax, and dt if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    if Hs0:
         ds["Hsig"]=ds["Hsig"]/Hs0
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    hs_levels = hs_levels or np.arange(0,ds["Hsig"].max().item(),ds["Hsig"].max().item()/100)
    hs_ticks = hs_ticks or np.around(np.arange(0,hs_levels.max()+hs_levels.max()/5,hs_levels.max()/5),decimals=3)
    # velocity clipping
    vel_clip_max=xr.apply_ufunc(np.abs,temp["Mfvx"]).quantile(vel_clip_max).item()
    temp["Mfvx"]=temp["Mfvx"].clip(min=-vel_clip_max,max=vel_clip_max)
    temp["Mfvy"]=temp["Mfvy"].clip(min=-vel_clip_max,max=vel_clip_max)
    vel_clip_min=xr.apply_ufunc(np.abs,temp["Mfvx"]).quantile(vel_clip_min).item()
    temp=temp.where(lambda x:(x["Mfvx"]**2+x["Mfvy"]**2)**0.5>=vel_clip_min,drop=True)
    # figure
    fig,ax=plt.subplots(figsize=(9.2,9.2*(ymax-ymin)/(xmax-xmin)))
    ax=[ax]
    ax[0].axis('equal')
    # plot hs
    hs=ds["Hsig"].plot(ax=ax[0],cmap=cmap,levels=hs_levels,add_colorbar=False,clip_on=True)
    fig.colorbar(hs,ticks=hs_ticks,label=[r'$\mathrm{ H_S}$ [m]',r'$\mathrm{ H_S\,/\,H_{S,0}}$'][Hs0 is not None],orientation='horizontal',cax=ax[0].inset_axes([0.05,1.05,0.85,0.05],transform=ax[0].transAxes),ticklocation='top')
    # plot land surface
    ax[0].fill_betweenx(ds.y,ds.isel(y=0).where(lambda x:x.Botlev<=0,drop=True).isel(x=0).x.item(),xmax,color="peachpuff",clip_on=True)
    # plot reef contour
    if plot_dep_dev:
        (-(ds.where(ds.Botlev-ds.isel(y=0).Botlev!=0,drop=True).Botlev)).plot.contourf(ax=ax[0],colors="k",add_colorbar=False)
    # plot quiver with Mfv
    quiv=temp.plot.quiver(ax=ax[0],x="x",y="y",u="Mfvx",v="Mfvy",scale=scale,add_guide=False)
    ax[0].quiverkey(quiv,0.95,1.02,vel_clip_max,f"{vel_clip_max:.2f} m/s")
    # plot depth contour
    if dep_levels is not None:
        x_dep_levels=[ds.isel(y=-1).where(lambda x:x["Botlev"]<=i,drop=True).isel(x=0).x.item() for i in dep_levels]
        depcont=(ds['Botlev']).plot.contour(levels=dep_levels,colors='grey',linewidth=3,linestyles="-",ax=ax[0])
        ax[0].clabel(depcont,fmt='-%.2f m',manual=[(i,ymin+(ymax-ymin)*0.9) for i in x_dep_levels],fontsize=14)
    # set axis properties
    [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]')) for i in ax]
    [(i.set_xlim([xmin,xmax]),i.set_ylim([ymin,ymax])) for i in ax]
    fig.tight_layout()
    plt.close()
    return fig

def Mfv_fig(ds,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=1,vel_clip_max=0.9,vel_clip_min=0.1,plot_dep_dev=False,dep_levels=None,mfvy_levels=None,mfvy_ticks=None,cmap="RdBu_r"):
    """
    Create 2D fig with mass flux velocities and y-component as colours.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'Mfvx' and Mfvy'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        vel_clip_max (float, optional): maximum x- and y-velocity clip (in % of x quantile). Default to 0.9.
        vel_clip_min (float, optional): minimum absolute velocity (in % of x quantile) to be plotted (otherwise nan). Default to 0.1
        plot_dep_dev (bool), optional: Condition plotting deviations from depth at the first cross-shore section. Default to False.        
        dep_levels(np array, optional): array with depth contour levels (in m) to be plotted. If None, no depth contours.
        mfvy_levels(np array, optional): array with V-currents contour levels (in m) to be plotted.
        mfvy_ticks(np array, optional): array with V-currents contour levels ticks (in m) to be plotted.
        cmap (str): matplotlib colour map. Default to "RdBu_r".
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, dx, ymin, ymax, and dy if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    mfvy_max=xr.apply_ufunc(np.abs,ds["Mfvy"]).max().item()
    mfvy_levels = mfvy_levels or np.arange(-mfvy_max,mfvy_max+mfvy_max/50,mfvy_max/50)
    mfvy_ticks = mfvy_ticks or np.around(np.arange(-mfvy_max,mfvy_max+mfvy_max/3,mfvy_max/3),decimals=3)
    # velocity clipping
    vel_clip_max=xr.apply_ufunc(np.abs,temp["Mfvx"]).quantile(vel_clip_max).item()
    temp["Mfvx"]=temp["Mfvx"].clip(min=-vel_clip_max,max=vel_clip_max)
    temp["Mfvy"]=temp["Mfvy"].clip(min=-vel_clip_max,max=vel_clip_max)
    vel_clip_min=xr.apply_ufunc(np.abs,temp["Mfvx"]).quantile(vel_clip_min).item()
    temp=temp.where(lambda x:(x["Mfvx"]**2+x["Mfvy"]**2)**0.5>=vel_clip_min,drop=True)
    # figure
    fig,ax=plt.subplots(figsize=(9.2,9.2*(ymax-ymin)/(xmax-xmin)))
    ax=[ax]
    # plot hs
    mfvy=ds["Mfvy"].plot(ax=ax[0],cmap=cmap,levels=mfvy_levels,add_colorbar=False,clip_on=True)
    fig.colorbar(mfvy,ticks=mfvy_ticks,label=r'$\mathrm{ V [m \cdot s^{-1}]}$ ',orientation='horizontal',cax=ax[0].inset_axes([0.05,1.05,0.85,0.05],transform=ax[0].transAxes),ticklocation='top')
    # plot land surface
    ax[0].fill_betweenx(ds.y,ds.isel(y=0).where(lambda x:x.Botlev<=0,drop=True).isel(x=0).x.item(),xmax,color="peachpuff",clip_on=True)
    # plot reef contour
    if plot_dep_dev:
        (-(ds.where(ds.Botlev-ds.isel(y=0).Botlev!=0,drop=True).Botlev)).plot.contourf(ax=ax[0],colors="k",add_colorbar=False)
    # plot quiver with Mfv
    quiv=temp.plot.quiver(ax=ax[0],x="x",y="y",u="Mfvx",v="Mfvy",scale=scale,add_guide=False)
    ax[0].quiverkey(quiv,0.95,1.02,vel_clip_max,f"{vel_clip_max:.2f} m/s")
    # plot depth contour
    if dep_levels is not None:
        x_dep_levels=[ds.isel(y=-1).where(lambda x:x["Botlev"]<=i,drop=True).isel(x=0).x.item() for i in dep_levels]
        depcont=(ds['Botlev']).plot.contour(levels=dep_levels,colors='grey',linewidth=3,linestyles="-",ax=ax[0])
        ax[0].clabel(depcont,fmt='-%.2f m',manual=[(i,ymin+(ymax-ymin)*0.9) for i in x_dep_levels],fontsize=14)
    # set axis properties
    [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]')) for i in ax]
    [(i.set_xlim([xmin,xmax]),i.set_ylim([ymin,ymax])) for i in ax]
    ax[0].axis('equal')
    fig.tight_layout()
    plt.close()
    return fig

def MV_fig(ds,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=1,vel_clip_max=0.9,vel_clip_min=0.1,plot_dep_dev=False,dep_levels=None,mveta_levels=None,mveta_ticks=None,cmap="RdBu_r"):
    """
    Create 2D fig with mass flux velocities and y-component as colours.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'MVksi' and MVeta'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        vel_clip_max (float, optional): maximum x- and y-velocity clip (in % of x quantile). Default to 0.9.
        vel_clip_min (float, optional): minimum absolute velocity (in % of x quantile) to be plotted (otherwise nan). Default to 0.1
        plot_dep_dev (bool), optional: Condition plotting deviations from depth at the first cross-shore section. Default to False.        
        dep_levels(np array, optional): array with depth contour levels (in m) to be plotted. If None, no depth contours.
        mveta_levels(np array, optional): array with V-currents contour levels (in m) to be plotted.
        mveta_ticks(np array, optional): array with V-currents contour levels ticks (in m) to be plotted.
        cmap (str): matplotlib colour map. Default to "RdBu_r".
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, dx, ymin, ymax, and dy if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    mveta_max=xr.apply_ufunc(np.abs,ds["MVeta"]).max().item()
    mveta_levels = mveta_levels or np.arange(-mveta_max,mveta_max+mveta_max/50,mveta_max/50)
    mveta_ticks = mveta_ticks or np.around(np.arange(-mveta_max,mveta_max+mveta_max/3,mveta_max/3),decimals=3)
    # velocity clipping
    vel_clip_max=xr.apply_ufunc(np.abs,temp["MVksi"]).quantile(vel_clip_max).item()
    temp["MVksi"]=temp["MVksi"].clip(min=-vel_clip_max,max=vel_clip_max)
    temp["MVeta"]=temp["MVeta"].clip(min=-vel_clip_max,max=vel_clip_max)
    vel_clip_min=xr.apply_ufunc(np.abs,temp["MVksi"]).quantile(vel_clip_min).item()
    temp=temp.where(lambda x:(x["MVksi"]**2+x["MVeta"]**2)**0.5>=vel_clip_min,drop=True)
    # figure
    fig,ax=plt.subplots(figsize=(9.2,9.2*(ymax-ymin)/(xmax-xmin)))
    ax=[ax]
    # plot hs
    mveta=ds["MVeta"].plot(ax=ax[0],cmap=cmap,levels=mveta_levels,add_colorbar=False,clip_on=True)
    fig.colorbar(mveta,ticks=mveta_ticks,label=r'$\mathrm{ V [m \cdot s^{-1}]}$ ',orientation='horizontal',cax=ax[0].inset_axes([0.05,1.05,0.85,0.05],transform=ax[0].transAxes),ticklocation='top')
    # plot land surface
    ax[0].fill_betweenx(ds.y,ds.isel(y=0).where(lambda x:x.Botlev<=0,drop=True).isel(x=0).x.item(),xmax,color="peachpuff",clip_on=True)
    # plot reef contour
    if plot_dep_dev:
        (-(ds.where(ds.Botlev-ds.isel(y=0).Botlev!=0,drop=True).Botlev)).plot.contourf(ax=ax[0],colors="k",add_colorbar=False)
    # plot quiver with MV
    quiv=temp.plot.quiver(ax=ax[0],x="x",y="y",u="MVksi",v="MVeta",scale=scale,add_guide=False)
    ax[0].quiverkey(quiv,0.95,1.02,vel_clip_max,f"{vel_clip_max:.2f} m/s")
    # plot depth contour
    if dep_levels is not None:
        x_dep_levels=[ds.isel(y=-1).where(lambda x:x["Botlev"]<=i,drop=True).isel(x=0).x.item() for i in dep_levels]
        depcont=(ds['Botlev']).plot.contour(levels=dep_levels,colors='grey',linewidth=3,linestyles="-",ax=ax[0])
        ax[0].clabel(depcont,fmt='-%.2f m',manual=[(i,ymin+(ymax-ymin)*0.9) for i in x_dep_levels],fontsize=14)
    # set axis properties
    [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]')) for i in ax]
    [(i.set_xlim([xmin,xmax]),i.set_ylim([ymin,ymax])) for i in ax]
    ax[0].axis('equal')
    fig.tight_layout()
    plt.close()
    return fig


def Hs_Theta_fig(ds,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=1,plot_dep_dev=False,dep_levels=None,hs_levels=None,hs_ticks=None,Hs0=None,cmap="jet"):
    """
    Create 2D fig with significant wave height and mass flux velocities.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'Hs', and 'Thetam'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        plot_dep_dev (bool), optional: Condition plotting deviations from depth at the first cross-shore section. Default to False.        
        dep_levels(np array, optional): array with depth contour levels (in m) to be plotted. If None, no depth contours.
        hs_levels(np array, optional): array with significant wave height contour levels (in m) to be plotted.
        hs_ticks(np array, optional): array with significant wave height contour levels ticks (in m) to be plotted.
        Hs0 (float, optional): deep water significant wave height (in m). If provided, Hs is normalised.
        cmap (str): hs matplotlib colour map. Default to "jet".
    """
    warnings.filterwarnings('ignore')
    ds["Thetax"]=np.cos(np.radians(ds["Thetam"]))
    ds["Thetay"]=np.sin(np.radians(ds["Thetam"]))    
    # Assign xmin, xmax, dx, ymin, ymax, dy, tmin, tmax, and dt if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    if Hs0:
         ds["Hsig"]=ds["Hsig"]/Hs0
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    hs_levels = hs_levels or np.arange(0,ds["Hsig"].max().item()+ds["Hsig"].max().item()/100,ds["Hsig"].max().item()/100)
    hs_ticks = hs_ticks or np.around(np.arange(0,hs_levels.max()+hs_levels.max()/5,hs_levels.max()/5),decimals=3)
    fig,ax=plt.subplots(figsize=(9.2,9.2*(ymax-ymin)/(xmax-xmin)))
    ax=[ax]
    ax[0].axis('equal')
    # plot hs
    hs=ds["Hsig"].plot(ax=ax[0],cmap=cmap,levels=hs_levels,add_colorbar=False,clip_on=True)
    fig.colorbar(hs,ticks=hs_ticks,label=[r'$\mathrm{ H_S}$ [m]',r'$\mathrm{ H_S\,/\,H_{S,0}}$'][Hs0 is not None],orientation='horizontal',cax=ax[0].inset_axes([0.05,1.05,0.85,0.05],transform=ax[0].transAxes),ticklocation='top')
    # plot land surface
    ax[0].fill_betweenx(ds.y,ds.isel(y=0).where(lambda x:x.Botlev<=0,drop=True).isel(x=0).x.item(),xmax,color="peachpuff",clip_on=True)
    # plot reef contour
    if plot_dep_dev:
        (-(ds.where(ds.Botlev-ds.isel(y=0).Botlev!=0,drop=True).Botlev)).plot.contourf(ax=ax[0],colors="k",add_colorbar=False)
    # plot quiver with Thetam
    quiv=temp.plot.quiver(ax=ax[0],x="x",y="y",u="Thetax",v="Thetay",scale=scale,add_guide=False)
    ax[0].quiverkey(quiv,0.95,1.02,1,"$\\theta_M$")
    # plot depth contour
    if dep_levels is not None:
        x_dep_levels=[ds.isel(y=-1).where(lambda x:x["Botlev"]<=i,drop=True).isel(x=0).x.item() for i in dep_levels]
        depcont=(ds['Botlev']).plot.contour(levels=dep_levels,colors='grey',linewidth=3,linestyles="-",ax=ax[0])
        ax[0].clabel(depcont,fmt='-%.2f m',manual=[(i,ymin+(ymax-ymin)*0.9) for i in x_dep_levels],fontsize=14)
    # set axis properties
    [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]')) for i in ax]
    [(i.set_xlim([xmin,xmax]),i.set_ylim([ymin,ymax])) for i in ax]
    fig.tight_layout()
    plt.close()
    return fig


def dt_fig(path_run,print="PRINT",Tp=None):
    """Plot time series of dt based on PRINT file

    Args:
        path_run (str): run path.
        print (str, optional): PRINT log file. Defaults to "PRINT".
        Tp (float,optional): If provided, show twin x axis at the top with time in Tp. Default to None.

    Returns:
        fig: matplotlib figure with time series of dt.
    """
    # import PRINT file
    file=[i.strip() for i in open(path_run+print,"r").readlines()]
    dt0=[float(i.split()[2]) for i in file if "COMPUTE" in i][0] # initial time step
    tinit=[float(i.split()[1][:2])*3600 + float(i.split()[1][2:4])*60 + float(i.split()[1][4:])  for i in file if "COMPUTE" in i][0] # initial time
    tend=[float(i.split()[4][:2])*3600 + float(i.split()[4][2:4])*60 + float(i.split()[4][4:])  for i in file if "COMPUTE" in i][0] # end time
    dt=[dt0]+[float(i.split("New time step:")[-1].split()[0]) for i in file if "New time step:" in i] # time steps array
    dt=np.repeat(np.array(dt),2)
    t=[[tinit]]+[[float(file[i].split("in sec:")[-1].split()[0]),float(file[i-2].split("in sec:")[-1].split()[0])] for i in range(2,len(file)) if "New time step:" in file[i-1]]+[[tend]] # time array
    t= np.array( [item for sublist in t for item in sublist])
    # figure
    fig,ax=plt.subplots(figsize=(10,5))
    ax.plot(t,dt,c="k",marker="o")
    xmin,xmax,dx=t[0],t[-1],200
    [(i.grid(),i.set_xlabel("t [s]"),i.set_ylabel("dt [s]"),i.set_xlim([xmax,xmax]),i.set_xticks(np.arange(xmin,xmax,dx))) for i in [ax]]
    [i.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) for i in [ax]]
    if Tp: # upper x-axis with time in Tp
        axup=ax.twiny()
        [(i.set_xticks(np.arange(xmin,xmax,dx)/Tp),i.set_xlim([xmin/Tp,xmax/Tp]),i.set_xlabel("t [Tp]")) for i in [axup]]    
    plt.close()
    return fig


def ener_enst_fig(path_run,ener_enstro="ener_enstro",Tp=None):
    """Plot time series of dt based on PRINT file

    Args:
        path_run (str): run path.
        ener_enstro (str, optional): file with ener_enstro. Defaults to "ener_enstro".
        Tp (float,optional): If provided, show twin x axis at the top with time in Tp. Default to None.

    Returns:
        fig: matplotlib figure with time series of energy and enstrophy.
    """
    # import energy enstrophy file (from swash*_further)
    ee=np.array([[float(i) for i in a.split()] for a in open(f'{path_run}/{ener_enstro}','r').readlines()[1:]])
    xmin,xmax,dx=0,np.around(ee[-1,0],-2)+100,200 #time axis properties
    # figure
    fig,ax=plt.subplots(nrows=2,figsize=(10,7),sharex=True)
    ax[0].plot(ee[:,0],ee[:,1],label="KE",c="k") #kinetic energy
    ax[0].plot(ee[:,0],ee[:,2],label="PE",c="r") #potential energy
    ax[1].plot(ee[:,0],ee[:,3],label="_nolegend_",c="k") #enstrophy
    [(i.grid(),i.set_xlim([xmin,xmax])) for i in ax]
    [i.set_ylabel(["$\mathrm{E\,[m^5 \cdot s^{-2}]}$","$\mathrm{Z\,[m^2 \cdot s^{-2}]}$"][j]) for j,i in enumerate(ax)]
    ax[-1].set_xlabel("$\mathrm{t\,[s]}$")
    [i.ticklabel_format(axis='y',style='sci',scilimits=(0,0)) for i in ax]
    [i.legend() for i in [ax[0]]]
    [ax[i].text(0.02,0.90,['a','b','c'][i],transform=ax[i].transAxes,ha='center',weight='bold') for i in range(len(ax))]
    ax[-1].set_xticks(np.arange(xmin,xmax,dx))
    if Tp: # upper x-axis with time in Tp
        axup=ax[0].twiny()
        [(i.set_xticks(np.arange(xmin,xmax,dx)/Tp),i.set_xlim([xmin/Tp,xmax/Tp]),i.set_xlabel("t [Tp]")) for i in [axup]]
    plt.close()
    return fig

def Hs_wm_fig(ds,ds_bot=None,dt=0.1,burst=30,Tp=2.83,Hs_wm=0.118):
    """Plot time series of modelled (total) Hs with prescribed time window.
    Note: later separate incoming and outgoing components.

    Args:
        ds (xr DataArray): Dataset with 'Watlev' (water levels at wave maker) (in m). If 'Vksi_k' (velocities at wave maker) (in m/s), calculate 
        ds_bot (xr DataArray): Dataset with 'Botlev' (bottom depth) (in m); only used for reflection calculation. Default to None (i.e., reflection is not calculated)
        dt (float, optional): time step (in s). Defaults to 0.1.
        burst (int, optional): burst length (in Tp). Defaults to 30.
        Tp (float, optional): peak wave period (in s). Defaults to 2.83.
        Hs_wm (float, optional): target significant wave height at the wave maker (in m). Defaults to 0.118.
    
    Returns:
        fig: matplotlib figure with time series of significant wave height at wave maker.        
    """
    warnings.filterwarnings('ignore')
    ds=ds.assign_coords(t=ds.t*dt) # assign time axis
    ds.t.attrs = {"standard_name": 'Tsec',"long_name":"Time in seconds from reference time","units":"s","axis":"t"}    
    xmin,xmax,dx=ds.t.min().item(),np.around(ds.t.max().item()),200 #time axis properties
    if "Vksi_k" in ds: # reflection
        ds["station_d"]=(("stations"),(ds_bot.Botlev.interp(x=ds.station_x,y=ds.station_y)).values) #calculate depth at station
        nV=len(ds.kc) #number of vertical layers
        ds=ds.assign_coords({"k":(1-(ds.kc/nV+1/(nV*2)))}) # create k coordinate (% height above bed; 0 for bed and 1 for water surface)
        ds["hab"]=(("stations","t","kc"),((ds.Watlev+ds.station_d)*ds.k).values,{"standard_name": "height above bed","units":"m"}) # from k = 0 (top) to k = nV-1 (near bed); velocities at mid grid
        ds=ds.isel(kc=-1,stations=0).drop("k") # bottom layer; this can be set as argument later
        etai,etar,ui,ur=lwt.etau_separation(t=ds.t.values,eta=ds.Watlev.values,u=ds.Vksi_k.values,h=(ds.Watlev.mean("t")+ds.station_d).item(),hab=ds.hab.mean("t").item(),fcutoff=0.5)
        ds["Watlev_i"]=("t",etai.real,{"standard_name": "incident water level","units":"m"})
        ds["Watlev_r"]=("t",etar.real,{"standard_name": "reflected water level","units":"m"})
        ds["Vksi_k_i"]=("t",ui.real,{"standard_name": "incident velocity","units":"m/s"})
        ds["Vksi_k_r"]=("t",ur.real,{"standard_name": "reflected velocity","units":"m/s"})
    # figure
    fig,ax=plt.subplots(figsize=(10,7))
    (4*ds.Watlev.rolling(t=int(Tp*burst/dt),center=True).std()).plot(ax=ax,c="k",label=f"Total Hs (window of {burst:.0f} Tps)")
    if "Vksi_k" in ds: # total
        (4*ds.Watlev_i.rolling(t=int(Tp*burst/dt),center=True).std()).plot(ax=ax,c="r",label=f"Hs_i (window of {burst:.0f} Tps)")
        (4*ds.Watlev_r.rolling(t=int(Tp*burst/dt),center=True).std()).plot(ax=ax,c="b",label=f"Hs_r (window of {burst:.0f} Tps)")
    ax.axhline(Hs_wm,c="k",ls="--",label="Target") #prescribed Hs
    [(i.set_xticks(np.arange(xmin,xmax,dx)),ax.set_ylabel("Hs [m] at wave maker"),i.grid(),i.legend(),i.set_xlim([xmin,xmax])) for i in [ax]]
    if Tp: # upper x-axis with time in Tp
        axup=ax.twiny()
        [(i.set_xticks(np.arange(xmin,xmax,dx)/Tp),i.set_xlim([xmin/Tp,xmax/Tp]),i.set_xlabel("t [Tp]")) for i in [axup]]
    plt.close()
    return fig

def unify_mkv(ds,mkmax=0.1):
    """Return dataset with merged 3D velocity profiles (Mkevx and Mkcvx --> Mkvx; Mkevy and Mkevy --> Mkvy).
    
    Args:
        ds (xr dataset): dataset with "Mkcvx", "Mkcvy", "Mkevx", "Mkevy" (in m/s).
        mkmax (float): top coordinate of water (in m) (command MKmax, see swash*_further).
        
    Returns:
        ds (xr dataset): dataset with "Mkvx", "Mkvy", and "z" variables added.
    """
    if "Mkcvy" not in ds: ds["Mkcvy"]=ds["Mkcvx"]*0
    if "Mkevy" not in ds: ds["Mkevy"]=ds["Mkcvx"]*0
    temp_zc=ds[["Mkcvx","Mkcvy"]]
    temp_zc=temp_zc.assign_coords({"kc":-(temp_zc["kc"]-0.5)/temp_zc["kc"].max()}).rename({"kc":"k","Mkcvx":"Mkvx","Mkcvy":"Mkvy"}) #transform kc into k; k is % of water depth, 0 for surface and -1 for bottom
    temp_ze=ds[["Mkevx","Mkevy"]]
    temp_ze=temp_ze.assign_coords({"ke":-(temp_ze["ke"])/temp_ze["ke"].max()}).rename({"ke":"k","Mkevx":"Mkvx","Mkevy":"Mkvy"}) #transform ke into k
    temp=xr.merge([temp_ze,temp_zc,ds[[i for i in ds.data_vars if i not in ["Mkcvx","Mkcvy","Mkevx","Mkevy"]]]])
    temp=temp.assign_coords({"z":mkmax+(temp.Botlev+mkmax)*temp.k}) #write z coordinate
    # rename attributes
    temp["Mkvx"].attrs["standard_name"],temp["Mkvy"].attrs["standard_name"]="Mkvx","Mkvy"
    temp["Mkvx"].attrs["long_name"],temp["Mkvy"].attrs["long_name"]="Mean layer-dependent u","Mean layer-dependent v"
    return temp

def trunc(x,decimals=3):
    """Truncate x to decimals."""
    x[x>0]=np.floor(x[x>0]*10**decimals)/10**decimals
    x[x<0]=np.ceil(x[x<0]*10**decimals)/10**decimals
    return x

def xvel_3D(ds,mkmax=0.1,xmin=None,xmax=None,dx=None,y=0,scale=2,quiver=True,line_vel=False,mfvx=False,mdavx=False,setup=False,contour=False,mkvx_levels=None,mkvx_ticks=None,cmap="RdBu_r"):
    """
    Plot cross-shore view with 2DV u-velocities.
    Args:
        ds (xr data structure): data structure with "Mkvx", "Mkvy" (in m/s), "Botlev" (in m), and optional parameters ( "Mdavx","Mdavy", "Mfvx", "Mfvy" (in m/s), and 'Setup' (m)).
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        dx (float, optional): x-grid size (m). If None, dx is calculated from ds.
        y (float, optional): y-position (in m). Default to 0.
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        quiver (bool, optional): if True, plot arrows with 2DV velocities. Default to True.
        line_vel (bool, optional): if True, plot vertical lines with 2DV velocities. Default to False.
        mfvx (bool, optional): if True, plot vertical dashed lines with mass-flux velocities. Default to False.
        mdavx (bool, optional): if True, plot vertical dashed lines with mean depth-averaged velocities. Default to False.
        setup (bool, optional): if True, plot horizontal dashed lines with setup. Default to False.
        contour (bool, optional): if True, plot colour map with "Mkvx". Default to False.
        mkvx_levels(np array, optional): array with "Mkvx" contour levels (in m) to be plotted.
        mkvx_ticks(np array, optional): array with "Mkvx" contour levels ticks (in m) to be plotted.
        cmap (str): matplotlib colour map. Default to "RdBu_r".

    Returns:
        fig: matplotlib figure with cross-shore view with 2DV u-velocities.
        """
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))))
    temp=temp.sel(y=y,method="nearest") # select y-position
    temp=unify_mkv(temp,mkmax=mkmax) # create condensed velocity profiles (interpolating Mkev and Mkcv; edges and centers)
    temp["Mkvz"]=temp["Mkvx"]*0 # create null z component
    # figure
    fig,ax=plt.subplots(figsize=(20,7))
    if contour:
        mkvx_max=xr.apply_ufunc(np.abs,temp["Mkvx"]).max().item()
        mkvx_levels = mkvx_levels or np.arange(-mkvx_max,mkvx_max+mkvx_max/50,mkvx_max/50)
        mkvx_ticks = mkvx_ticks or trunc(np.arange(-mkvx_max,mkvx_max+mkvx_max/2,mkvx_max),decimals=3)
        mkvx=temp.Mkvx.plot(x='x',y='z',cmap=cmap,levels=mkvx_levels,add_colorbar=False)
        fig.colorbar(mkvx,ticks=mkvx_ticks,label=r'$\mathrm{ u\, [m \cdot s^{-1}]}$ ',orientation='horizontal',cax=ax.inset_axes([0.45,0.05,0.35,0.03],transform=ax.transAxes),ticklocation='top') # depth profile from left to right 
    xarray=np.arange(xmin,xmax,dx) # create area of plot
    for x in xarray:
        inst=temp.sel(x=x,method="nearest") #select data for given x
        if line_vel: ax.plot(inst.x+4*scale*inst.Mkvx,inst.z,c="k",marker="o") #draw 3D vel profile; I don't get why I need 5*scale
        if quiver:
            quiv=inst.plot.quiver(x="x",y="z",u="Mkvx",v="Mkvz",ax=ax,pivot="tail",scale=scale,add_guide=False,width=0.002,headwidth=5)
            ax.quiverkey(quiv,0.95,1.02,0.1,f"{0.1:.2f} m/s")
        ax.plot([inst.x]*2,[(-inst.Botlev),mkmax],c="k",ls="--") #draw ref line
        if mfvx: ax.plot([inst.x+inst.Mfvx*5*scale]*2,[(-inst.Botlev),mkmax],c="r",ls="--") #draw Mfvx
        if mdavx: ax.plot([inst.x+inst.Mdavx*5*scale]*2,[(-inst.Botlev),mkmax],c="b",ls="--") #draw Mdavx
    (-ds.Botlev).sel(x=slice(xmin,xmax)).sel(y=y,method="nearest").plot(c="k",label="Bottom") # plot depth profile
    if setup: ds.sel(x=slice(xmin,xmax)).sel(y=y,method="nearest").where(lambda ds:ds["Setup"]>-ds["Botlev"]).Setup.plot(c="r",label="Setup") # plot setup profile
    [(i.set_xlabel("x [m]"),i.set_ylabel("z [m]")) for i in [ax]] #i.set_xlim([xmin,xmax])
    [(i.grid(),i.axhline(0,c="grey",ls="-"),i.legend()) for i in [ax]]
    ax.text(0.5,1.1,"2DV velocity profile",transform=ax.transAxes,fontsize=14,ha="center")
    plt.close()
    return fig

def yvel_3D(ds,mkmax=0.1,ymin=None,ymax=None,dy=None,x=0,scale=2,quiver=True,line_vel=False,mfvy=False,mdavy=False,setup=False,contour=False,mkvy_levels=None,mkvy_ticks=None,cmap="RdBu_r"):
    """
    Plot alongshore view with 2DV v-velocities.
    Args:
        ds (xr data structure): data structure with "Mkvx", "Mkvy", "Mdavx","Mdavy", "Mfvx", "Mfvy" (in m/s), and "Botlev" (in m).
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        dy (float, optional): y-grid size (m). If None, dy is calculated from ds.
        x (float, optional): x-position (in m). Default to 0.
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        quiver (bool, optional): if True, plot arrows with 2DV velocities. Default to True.
        line_vel (bool, optional): if True, plot vertical lines with 2DV velocities (on top of arrows). Default to False.
        mfvy (bool, optional): if True, plot vertical dashed lines with mass-flux velocities. Default to False.
        mdavy (bool, optional): if True, plot vertical dashed lines with mean depth-averaged velocities. Default to False.
        setup (bool, optional): if True, plot horizontal dashed lines with setup. Default to False.
        contour (bool, optional): if True, plot colour map with "Mkvy". Default to False.
        mkvy_levels(np array, optional): array with "Mkvy" contour levels (in m) to be plotted.
        mkvy_ticks(np array, optional): array with "Mkvy" contour levels ticks (in m) to be plotted.
        cmap (str): matplotlib colour map. Default to "RdBu_r".

    Returns:
        fig: matplotlib figure with cross-shore view with 2DV u-velocities.
        """
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()
    temp=ds.sel(y=slice(ymin,ymax,math.floor(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    temp=temp.sel(x=x,method="nearest") # select x-position
    temp=unify_mkv(temp,mkmax=mkmax) # condenses velocity profiles (interpolating Mkev and Mkcv; edges and centers)
    temp["Mkvz"]=temp["Mkvx"]*0 # create null z component
    # figure
    fig,ax=plt.subplots(figsize=(20,7))
    if contour:
        mkvy_max=xr.apply_ufunc(np.abs,temp["Mkvy"]).max().item()
        mkvy_levels = mkvy_levels or np.arange(-mkvy_max,mkvy_max+mkvy_max/50,mkvy_max/50)
        mkvy_ticks = mkvy_ticks or trunc(np.arange(-mkvy_max,mkvy_max+mkvy_max/2,mkvy_max),decimals=3)
        mkvy=temp.Mkvy.plot(x='y',y='z',cmap=cmap,levels=mkvy_levels,add_colorbar=False)
        fig.colorbar(mkvy,ticks=mkvy_ticks,label=r'$\mathrm{ v\, [m \cdot s^{-1}]}$ ',orientation='horizontal',cax=ax.inset_axes([0.45,0.05,0.35,0.03],transform=ax.transAxes),ticklocation='top') # depth profile from left to right 
    yarray=np.arange(ymin,ymax,dy) # create area of plot
    for y in yarray:
        inst=temp.sel(y=y,method="nearest") #select data for given y
        if line_vel: ax.plot(inst.y+5*scale*inst.Mkvy,inst.z,c="k",marker="o") #draw 3D vel profile; I don't get why I need 5*scale
        if quiver:
            quiv=inst.plot.quiver(x="y",y="z",u="Mkvy",v="Mkvz",ax=ax,pivot="tail",scale=scale,add_guide=False,width=0.002,headwidth=5)
            ax.quiverkey(quiv,0.95,1.02,0.1,f"{0.1:.2f} m/s")
        ax.plot([inst.y]*2,[(-inst.Botlev),mkmax],c="k",ls="--") #draw ref line
        if mfvy: ax.plot([inst.y+inst.Mfvy*5*scale]*2,[(-inst.Botlev),mkmax],c="r",ls="--") #draw Mfvy
        if mdavy: ax.plot([inst.y+inst.Mdavy*5*scale]*2,[(-inst.Botlev),mkmax],c="b",ls="--") #draw Mdavx
    (-ds.Botlev).sel(y=slice(ymin,ymax)).sel(x=x,method="nearest").plot(c="k",label="Bottom") # plot depth profile
    if setup: ds.sel(y=slice(ymin,ymax)).sel(x=x,method="nearest").where(lambda ds:ds["Setup"]>-ds["Botlev"]).Setup.plot(c="r",label="Setup") # plot setup  profile
    ax.axhline(0,c="grey",ls="-")
    ax.axvline(0,c="grey",ls="-")
    ax.grid()
    ax.text(0.5,1.1,"2DV velocity profile",transform=ax.transAxes,fontsize=14,ha="center")
    plt.close()
    return fig


def ibp_nan(ds,exc_value=10):
    """
    Plot alongshore view with 2DV v-velocities.
    Args:
        ds (xr data structure): data structure with "Ibp" (in m).
        exc_value (float, optional): Exception value for runup (abs values greater than exc_value are excluded; in m). Default to 10.

        Returns:
            fig: matplotlib figure with % of nan values.
        """
    fig,ax=plt.subplots(figsize=(10,7))
    temp=ds.where(lambda x:np.abs(x["Ibp"])>exc_value)
    temp["Ibp_nan"]=(xr.apply_ufunc(np.isnan,temp["Ibp"])).sum(dim="t")/temp["Ibp"].isel(y=0).size*100
    temp["Ibp_nan"].plot(ax=ax,c="k",clip_on=True)
    ax.set_ylabel("% of nan values in Ibp")
    ax.grid()
    ax.set_ylim([-1,101])
    plt.close()
    return fig

def runup_fig(R2,Ibp,R2_S06,max_z,exc_value=10,ymin=None,ymax=None,dy=None,Lr=None):
    """Plot along-section of runup properties.

    Args:
        R2 (xr DataArray): xr dataarray with 2% runup (in m).
        Ibp (xr DataArray): xr dataarray with Ibp (in m).
        R2_S06 (float): R2 according to Stockdon (2006).
        max_z (float): maximum bottom level (in m)
        exc_value (float, optional): Exception value for runup (abs values greater than exc_value are excluded; in m). Default to 10.
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        dy (float, optional): y-grid size (m). If None, dy is calculated from ds.
        Lr (float, optional): reef lenght (m). If provided, draw vertical lines at reef position.

    Return:
        figure with along-section of runup properties.
    """
    ymin = ymin or R2.y.min().item()
    ymax = ymax or R2.y.max().item()
    dy = dy or (R2.y.isel(y=1)-R2.y.isel(y=0)).item()
    R2=R2.sel(y=slice(ymin,ymax,math.floor(dy/((R2.y.isel(y=1)-R2.y.isel(y=0)).values.item()))))    
    Ibp=Ibp.sel(y=slice(ymin,ymax,math.floor(dy/((Ibp.y.isel(y=1)-Ibp.y.isel(y=0)).values.item()))))    
    if exc_value:
        bad_data=(xr.apply_ufunc(np.isnan,Ibp.where(lambda x:np.abs(x)>exc_value)).sum()/Ibp.size*100).values.item()
        if bad_data!=100:
            print(f'{100-bad_data:.1f} % of runup data is corrupted')
        Ibp=Ibp.where(lambda x:np.abs(x)<exc_value,drop=True)    
    fig,ax=plt.subplots(figsize=(10,7))
    R2.plot(c="k",label="R2")
    (Ibp.mean(dim="t")+2*Ibp.std(dim="t")).plot(c="k",label="R2_spec",ls="--")
    ax.axhline(y=R2_S06,c="k",label="R2_S06",ls="-.")
    Ibp.max(dim="t").plot(c="g",label="Max BWL")
    Ibp.min(dim="t").plot(c="g",ls="--",label="Min BWL")
    Ibp.mean(dim="t").plot(c="r",label="Setup")
    ax.axhline(y=max_z,c="b",label="Max z")
    if Lr:
        [ax.axvline(x=i*Lr/2,c="k",ls="--") for i in [1,-1]]
    [(i.grid(),i.set_ylabel("z [m]"),i.legend(ncol=7),i.set_title(""),i.set_xlim([ymin,ymax])) for i in [ax]]
    plt.close()
    return fig

def cprof_fig(ds1,ds2,ds3,y=[0,0,-10],c=["r","b","k"],height=200):
    """Return figure with cross-shore view of Hs, Setup, U, and depth profiles for 3 runs (ds1-3).

    Args:
        ds1 (xr Dataset): xr Dataset 1 with 'Hsig', 'Setup', 'MVksi' and 'Botlev'.
        ds2 (xr Dataset): xr Dataset 2 with 'Hsig', 'Setup', 'MVksi' and 'Botlev'.
        ds3 (xr Dataset): xr Dataset 3 with 'Hsig', 'Setup', 'MVksi' and 'Botlev'.
        y (list, optional): list with y-values where to plot cross-section for ds1, ds2, and ds3. Defaults to [0,0,-10].
        c (list, optional): list with color strings of ds1, ds2, and ds3, respectively. Defaults to ["r","b","k"].
        height (int, optional): _description_. Defaults to 200.

    Returns:
        _type_: _description_
    """
    import holoviews as hv
    return hv.Layout(ds1.sel(y=y[0]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Hsig.hvplot(c=c[0]).opts(gridstyle={'grid_line_color': 'grey'}, show_grid=True,height=height,title='',xlabel='',xticks='',ylabel="Hs [m]")
    *\
    ds2.sel(y=y[1]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Hsig.hvplot(c=c[1]) *\
    ds3.sel(y=y[2]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Hsig.hvplot(c=c[2]) +\
    ds1.sel(y=y[0]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Setup.hvplot(c=c[0]).opts(gridstyle={'grid_line_color': 'grey'}, show_grid=True,height=height,title='',xlabel='',ylabel="$\eta$ [m]") *\
    ds2.sel(y=y[1]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Setup.hvplot(c=c[1]) *\
    ds3.sel(y=y[2]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Setup.hvplot(c=c[2])+\
    ds1.sel(y=y[0]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).MVksi.hvplot(c=c[0]).opts(gridstyle={'grid_line_color': 'grey'}, show_grid=True,height=height,title='',xlabel='',ylabel="U [m/s]") *\
    ds2.sel(y=y[1]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).MVksi.hvplot(c=c[1]) *\
    ds3.sel(y=y[2]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).MVksi.hvplot(c=c[2])+\
    (-ds1.sel(y=y[0]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Botlev).hvplot(c=c[0]).opts(gridstyle={'grid_line_color': 'grey'}, show_grid=True,height=height,title='',ylabel="d [m]") *\
    (-ds2.sel(y=y[1]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Botlev).hvplot(c=c[1]) *\
    (-ds3.sel(y=y[2]).where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True).Botlev).hvplot(c=c[2])).cols(1)
    
def Mamfv_fig(ds,xmin=None,xmax=None,dx=None,ymin=None,ymax=None,dy=None,scale=1,vel_clip_max=0.9,vel_clip_min=0.1,plot_dep_dev=False,dep_levels=None,mamfvy_levels=None,mamfvy_ticks=None,cmap="RdBu_r"):
    """
    Create 2D fig with moving average of mass flux velocities and y-component as colours.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'Mamfvx' and Mamfvy'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        scale (float, optional): quiver scale (larger values result in smaller arrows). Default to 1.
        vel_clip_max (float, optional): maximum x- and y-velocity clip (in % of x quantile). Default to 0.9.
        vel_clip_min (float, optional): minimum absolute velocity (in % of x quantile) to be plotted (otherwise nan). Default to 0.1
        plot_dep_dev (bool), optional: Condition plotting deviations from depth at the first cross-shore section. Default to False.        
        dep_levels(np array, optional): array with depth contour levels (in m) to be plotted. If None, no depth contours.
        mamfvy_levels(np array, optional): array with V-currents contour levels (in m) to be plotted.
        mamfvy_ticks(np array, optional): array with V-currents contour levels ticks (in m) to be plotted.
        cmap (str): matplotlib colour map. Default to "RdBu_r".
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, dx, ymin, ymax, and dy if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    dx = dx or (ds.x.isel(x=1)-ds.x.isel(x=0)).item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    dy = dy or (ds.y.isel(y=1)-ds.y.isel(y=0)).item()    
    temp=ds.sel(x=slice(xmin,xmax,math.ceil(dx/((ds.x.isel(x=1)-ds.x.isel(x=0)).values.item()))),
                y=slice(ymin,ymax,math.ceil(dy/((ds.y.isel(y=1)-ds.y.isel(y=0)).values.item()))))
    mamfvy_max=xr.apply_ufunc(np.abs,ds["Mamfvy"]).max().item()
    if mamfvy_levels is None: mamfvy_levels = np.arange(-mamfvy_max,mamfvy_max+mamfvy_max/50,mamfvy_max/50)
    if mamfvy_ticks is None: mamfvy_ticks = np.around(np.arange(-mamfvy_max,mamfvy_max+mamfvy_max/3,mamfvy_max/3),decimals=3)
    # velocity clipping
    vel_clip_max=xr.apply_ufunc(np.abs,temp["Mamfvx"]).quantile(vel_clip_max).item()
    temp["Mamfvx"]=temp["Mamfvx"].clip(min=-vel_clip_max,max=vel_clip_max)
    temp["Mamfvy"]=temp["Mamfvy"].clip(min=-vel_clip_max,max=vel_clip_max)
    vel_clip_min=xr.apply_ufunc(np.abs,temp["Mamfvx"]).quantile(vel_clip_min).item()
    temp=temp.where(lambda x:(x["Mamfvx"]**2+x["Mamfvy"]**2)**0.5>=vel_clip_min,drop=True)
    # figure
    fig,ax=plt.subplots(figsize=(9.2,9.2*(ymax-ymin)/(xmax-xmin)))
    ax=[ax]
    # plot Mamfv
    mamfvy=ds["Mamfvy"].plot(ax=ax[0],cmap=cmap,levels=mamfvy_levels,add_colorbar=False,clip_on=True)
    fig.colorbar(mamfvy,ticks=mamfvy_ticks,label=r'$\mathrm{ V [m \cdot s^{-1}]}$ ',orientation='horizontal',cax=ax[0].inset_axes([0.05,1.05,0.85,0.05],transform=ax[0].transAxes),ticklocation='top')
    # plot land surface
    ax[0].fill_betweenx(ds.y,ds.isel(y=0).where(lambda x:x.Botlev<=0,drop=True).isel(x=0).x.item(),xmax,color="peachpuff",clip_on=True)
    # plot reef contour
    if plot_dep_dev:
        (-(ds.where(ds.Botlev-ds.isel(y=0).Botlev!=0,drop=True).Botlev)).plot.contourf(ax=ax[0],colors="k",add_colorbar=False)
    # plot quiver with Mamfv
    quiv=temp.plot.quiver(ax=ax[0],x="x",y="y",u="Mamfvx",v="Mamfvy",scale=scale,add_guide=False)
    ax[0].quiverkey(quiv,0.95,1.02,vel_clip_max,f"{vel_clip_max:.2f} m/s")
    # plot depth contour
    if dep_levels is not None:
        x_dep_levels=[ds.isel(y=-1).where(lambda x:x["Botlev"]<=i,drop=True).isel(x=0).x.item() for i in dep_levels]
        depcont=(ds['Botlev']).plot.contour(levels=dep_levels,colors='grey',linewidth=3,linestyles="-",ax=ax[0])
        ax[0].clabel(depcont,fmt='-%.2f m',manual=[(i,ymin+(ymax-ymin)*0.9) for i in x_dep_levels],fontsize=14)
    # set axis properties
    [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]')) for i in ax]
    [(i.set_xlim([xmin,xmax]),i.set_ylim([ymin,ymax])) for i in ax]
    ax[0].axis('equal')
    fig.tight_layout()
    plt.close()
    return fig

def crossp_fig(ds,y_reef,y_exposed,xmin=None,xmax=None):
    """
    Create 2D fig with moving average of mass flux velocities and y-component as colours.
    Args:
        ds (xr data structure): Single data structure with 'x', 'Botlev', 'Hsig', 'Setup', 'Mdavx', and 'Mdavy'.
        y1 (float): position 1 of y where to plot cross-section.
        y2 (float): position 1 of y where to plot cross-section.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, dx, ymin, ymax, and dy if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    temp=ds.sel(x=slice(xmin,xmax))
    # figure
    fig,ax=plt.subplots(figsize=(10,15),nrows=5,sharex=True)
    [temp.sel(y=[y_reef,y_exposed][j],method="nearest").where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True)["Hsig"].interpolate_na('x').plot(ax=ax[0],c=["r","k"][j],label=["Reef","Exposed"][j]) for j in range(2)]
    [temp.sel(y=[y_reef,y_exposed][j],method="nearest").where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True)["Setup"].interpolate_na('x').plot(ax=ax[1],c=["r","k"][j],label=["Reef","Exposed"][j]) for j in range(2)]
    [temp.sel(y=[y_reef,y_exposed][j],method="nearest").where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True)["Mdavx"].interpolate_na('x').plot(ax=ax[2],c=["r","k"][j],label=["Reef","Exposed"][j]) for j in range(2)]    
    [temp.sel(y=[y_reef,y_exposed][j],method="nearest").where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True)["Mdavy"].interpolate_na('x').plot(ax=ax[3],c=["r","k"][j],label=["Reef","Exposed"][j]) for j in range(2)]    
    [(-temp.sel(y=[y_reef,y_exposed][j],method="nearest").where(lambda ds:ds["Setup"]>=-ds["Botlev"],drop=True)["Botlev"]).plot(ax=ax[-1],c=["r","k"][j],label=["Reef","Exposed"][j]) for j in range(2)]        
    [ax[i].set_ylabel(["$\mathrm{H_S\,\,[m]}$","$\mathrm{\eta\,\,[m]}$","$\mathrm{\\langle U_{dep}\\rangle\,\,[m \cdot s^{-1}]}$","$\mathrm{\\langle V_{dep}\\rangle\,\,[m \cdot s^{-1}]}$","$\mathrm{d\,\,[m]}$"][i]) for i in range(len(ax))]
    [i.set_xlabel("") for i in ax[:-1]]
    ax[-1].set_xlabel("x [m]")
    [(i.grid(),i.set_title(""),i.set_xlim([xmin,xmax])) for i in ax]
    ax[0].legend(ncol=3,bbox_to_anchor=(0.5,1.06),loc="lower center")
    [j.text(0,1.02,string.ascii_lowercase[i],fontweight="bold",transform=j.transAxes) for i,j in enumerate(ax.flatten())]
    plt.close()
    return fig

if __name__ == '__main__':
    pass

