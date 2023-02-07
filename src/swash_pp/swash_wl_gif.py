import warnings, math, os
import matplotlib.pyplot as plt
import imageio
import xarray as xr
import matplotlib.colors as colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

def create_2D_gif(gif_name,ds,xmin=None,xmax=None,tmin=None,tmax=None,dt=None,zmin=None,zmax=None,depmin=None,acc_factor=1,gif_path="",axis_off=False,dpi=100):
    """
    Create 2D gif with water level animation
    Args:
        gif_name (str): Gif name
        ds (xr data structure): Single data structure with 'x', 'Botlev' and 'Watlev'.
        xmin (float): minimum x-position (m). If None, ds.x.min().
        xmax (float): maximum x-position (m). If None, ds.x.max().
        tmin (float): minimum time (s). If None, ds.t.min().
        tmax (float): maximum time (s). If None, ds.t.max().
        dt (float): time step (s). If None, (ds.t.isel(t=1)-ds.t.isel(t=0)).
        zmin (float): minimum z-position (m). If None, (-ds.Botlev).min().
        zmax (float): maximum z-position (m). If None, ds.Watlev.max().
        depmin (float, optional): minimum threshold depth for runup calculation (m). Default to None.
        acc_factor (float, optional): Ratio of frame duration by dt. Default to 1.
        gif_path (str, optional): path where to save gif. Default to "" (local path).
        axis_off (bool, optional): boolean for plotting axis. Default to False (i.e., axis is plotted).
        dpi (int): Image dpi. Default to 100.
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, tmin, tmax, dt, zmin, and xmax if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    tmin = tmin or ds.t.min().item()
    tmax = tmax or ds.t.max().item()
    dt = dt or (ds.t.isel(t=1)-ds.t.isel(t=0)).item()
    temp=ds.sel(x=slice(xmin,xmax),t=slice(tmin,tmax,math.ceil(dt/((ds.t.isel(t=1)-ds.t.isel(t=0)).values.item()))))
    zmin = zmin or (-temp.Botlev).min().item()
    zmax = zmax or max(temp.Watlev.max().item(),(-temp.Botlev).max().item())
    if not os.path.exists(gif_name): os.mkdir(f'{gif_path}{gif_name}')
    d=temp.Botlev # bottom level
    wl=temp.Watlev # water level
    if depmin: # for runup calaculation
        temp["Dep"]=temp["Botlev"]+temp["Watlev"]
        xerror=-1000 #negative value for runup calculation
        xrunup=xr.where(temp["Dep"]>depmin,temp["x"],xerror).max(dim="x") # to locate waterline position
        ibp=temp.sel(x=xrunup,t=xrunup.t).Watlev # create instantaneous beach position array (or water level)
        # calculate runup range
        vmin,vmax = ibp.min().item(),ibp.max().item()
        if vmin>0: vmin=-(vmax-vmin)*.2
    def fig_rp(t=0):
        _,ax=plt.subplots(constrained_layout=True,figsize=(10,3))
        ax=[ax]
        # plot land and water level surfaces
        ax[0].fill_between(d.x,zmin+0*d.values.squeeze(),-d.values.squeeze(),color="peachpuff") # contour beach
        ax[0].fill_between(d.x,-d.values.squeeze(),wl.isel(t=t).values.squeeze(),color="skyblue",alpha=0.5)
        if depmin: # runup line and scatter
            aax=ax[0].inset_axes(
                [0.7,0.2,0.2,0.2], transform=ax[0].transAxes)
            if axis_off:
                [(i.set_xticklabels(""),i.set_yticklabels(""),i.axis([ibp.t.isel(t=0),ibp.t.isel(t=-1),vmin,vmax]),i.set_xlabel("t"),i.set_ylabel("z (SWL)"),i.grid("on"),i.axhline(y=0,c="k",ls="--",lw=2)) for i in [aax]]
            else:
                [(i.axis([ibp.t.isel(t=0),ibp.t.isel(t=-1),vmin,vmax]),i.set_xlabel("t [s]"),i.set_ylabel("z [m]"),i.grid("on"),i.axhline(y=0,c="k",ls="--",lw=2)) for i in [aax]]
            aax.plot([ibp.isel(t=t).t.item()],[ibp.isel(t=t).item()],ls="",marker="o",c="r")
            aax.plot(ibp.isel(t=slice(0,t+1)).t.values,ibp.isel(t=slice(0,t+1)).values,ls="-",marker="",c="r")
            ax[0].scatter(xrunup.isel(t=t),ibp.isel(t=t),c="r",s=50)
            ax[0].plot([xrunup.isel(t=t).item()]*2,[zmin,ibp.isel(t=t).item()],c="r",ls="--") # at SWL
        # set axis properties
        ax[0].set_title("")
        if axis_off:
            ax[0].axis('off')
        else:
            [(i.set_xlabel('X [m]'),i.set_ylabel('$\zeta,-d $ [m]')) for i in ax]
            ax[0].plot([(-d).where(lambda x:x>=0,drop=True).isel(x=0).x.item(),xmax],[0,0],c="k",ls="--") # at SWL
            ax[0].text(0.5,1.1,'t = {:.1f} s'.format(temp.isel(t=t).t.item()),transform=ax[0].transAxes,ha='center',va='bottom') # print time
        ax[0].axhline(c="k",ls="--") # at SWL
        ax[0].text(xmax,0,"SWL",ha="right")
        [(i.set_xlim([xmin,xmax]),i.set_ylim([zmin,zmax])) for i in ax]
        plt.savefig(f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png',dpi=dpi)
        plt.close() 

    print(f"Creating figures - total of {len(temp.t)} time steps")
    for t in range(len(temp.t)):
        if t%50==0: print(f"{t+1}/{len(temp.t)}")
        fig_rp(t=t)
        
    print(f"Creating gif")
    frames_gif(gf=len(temp.t),dt=dt,gif_name=gif_name,gif_path=gif_path,acc_factor=acc_factor)
    print(f"Deleting figures")
    delete_fig(gf=len(temp.t),gif_name=gif_name,gif_path=gif_path)


def create_3D_gif(gif_name,ds,xmin=None,xmax=None,ymin=None,ymax=None,tmin=None,tmax=None,dt=None,zmin=None,zmax=None,
                  aspect_ratio=None,
                  vmin=None,vmax=None,sc_wlv=1,elev=50,elev_light=1,azim=-135,azim_light=-155,depmin=None,
                  acc_factor=1,gif_path="",axis_off=False,dpi=100,combine_only=False,tpar=-1):
    """
    Create 3D gif with water level animation
    Args:
        gif_name (str): Gif name
        ds (xr data structure): Single data structure with 'x', 'Botlev' and 'Watlev'.
        xmin (float, optional): minimum x-position (m). If None, ds.x.min().
        xmax (float, optional): maximum x-position (m). If None, ds.x.max().
        ymin (float, optional): minimum y-position (m). If None, ds.y.min().
        ymax (float, optional): maximum y-position (m). If None, ds.y.max().
        tmin (float, optional): minimum time (s). If None, ds.t.min().
        tmax (float, optional): maximum time (s). If None, ds.t.max().
        dt (float, optional): time step (s). If None, (ds.t.isel(t=1)-ds.t.isel(t=0)).
        zmin (float, optional): minimum z-position (m). If None, (-ds.Botlev).min().
        zmax (float, optional): maximum z-position (m). If None, (-ds.Botlev).max().
        aspect_ratio (list, optional): [x,y,z] aspect ratio. If None, [15,7,5].
        vmin (float, optional): minimum water level (m). If None, ds.Watlev.min().
        vmax (float, optional): maximum water level (m). If None, ds.Watlev.max().
        sc_wlv (float, optional): scale of water level in 3D plot. Default to 1.   
        elev (float, optional): elevation view (m). Default to 50.
        elev_light (float, optional): elevation view light (m). Default to 50.
        azim (float, optional): azimuth view (deg). Default to -135.
        azim_light (float, optional): azimuth view light (deg). Default to -155.
        depmin (float, optional): minimum threshold depth for runup calculation (m). Default to None.
        acc_factor (float, optional): Ratio of frame duration by dt. Default to 1.
        gif_path (str, optional): path where to save gif. Default to "" (local path).
        axis_off (bool, optional): boolean for plotting axis. Default to False (i.e., axis is plotted).
        dpi (int): Image dpi. Default to 100.
        combine_only (bool): Condition for not generating frames but only combining them. Default to False.
        tpar (int): If different than -1, only generate frame for index tpar (for parallel computation). Default to -1.
    
        For parallel computing, the argument tpar can be used across a parallelized time loop. For example:
    '''
    # define the parallel function
    import multiprocessing
    def process_elements_in_parallel(time):
        with multiprocessing.Pool(processes=4) as pool: # 4 is the number of CPUs to be used
            pool.starmap(wlg.create_3D_gif, [(gif_name,ds,xmin,xmax,ymin,ymax,tmin,tmax,dt,zmin,zmax,
                    aspect_ratio,
                    vmin,vmax,sc_wlv,elev,elev-light,azim,azim_light,depmin,
                    acc_factor,gif_path,axis_off,dpi,False,t) for t in time]) #note that here this function needs all arguments to be explicitly stated    
    # with tmin, tmax, and dt, define temp xr dataset and get the time length
    t_len =len(ds.sel(t=slice(tmin,tmax,math.ceil(dt/((ds.t.isel(t=1)-ds.t.isel(t=0)).values.item())))).t)
    # generate frames in parallel
    process_elements_in_parallel(range(t_len)) 
    # combine frames and create 3D gif
    wlg.create_3D_gif(gif_name=gif_name,ds=ds,tmin=tmin,tmax=tmax,dt=dt,acc_factor=acc_factor,combine_only=True) 
    '''
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, ymin, ymax, tmin, tmax, dt, zmin, zmax, aspect_ratio, vmin, vmax, and sc_wlv, if not defined
    xmin = xmin or ds.x.min().item()
    xmax = xmax or ds.x.max().item()
    ymin = ymin or ds.y.min().item()
    ymax = ymax or ds.y.max().item()    
    tmin = tmin or ds.t.min().item()
    tmax = tmax or ds.t.max().item()
    dt = dt or (ds.t.isel(t=1)-ds.t.isel(t=0)).item()
    # create subset of xarray dataset
    temp=ds.sel(x=slice(xmin,xmax),y=slice(ymin,ymax),t=slice(tmin,tmax,math.ceil(dt/((ds.t.isel(t=1)-ds.t.isel(t=0)).values.item()))))
    zmin = zmin or (-temp.Botlev).min().item()
    zmax = zmax or (-temp.Botlev).max().item()
    aspect_ratio = aspect_ratio or [15,7,5]
    vmin = vmin or temp.Watlev.min().item()
    vmax = vmax or temp.Watlev.max().item()
    if not os.path.exists(gif_name): os.mkdir(f'{gif_path}{gif_name}')
    if depmin: # for runup calaculation
        temp["Dep"]=temp["Botlev"]+temp["Watlev"]
        xerror=-1000 #negative value for runup calculation
        xrunup=xr.where(temp["Dep"]>depmin,temp["x"],xerror).max(dim="x") # to locate waterline position
        temp["Ibp"]=temp.sel(x=xrunup,t=xrunup.t).Watlev # create instantaneous beach position array (or water level)
        # calculate runup range
        ibpmin,ibpmax=temp["Ibp"].min().item(),temp["Ibp"].max().item()
        if ibpmin>0: ibpmin=-(ibpmax-ibpmin)*.2
    #zorder: seabed,water surf,dep centre,runup
    zorder=[30,40,5,60]
    def wlv_1d(ax,t=0): # upper 2D plot
        # plot land and water level surfaces
        ax.fill_between(temp.x,zmin+0*(-temp.Botlev.isel(y=len(temp.y)//2)).values.squeeze(),(-temp.Botlev.isel(y=len(temp.y)//2)).values.squeeze(),color="peachpuff") # contour beach
        ax.fill_between(temp.x,(-temp.Botlev.isel(y=len(temp.y)//2)).values.squeeze(),temp.Watlev.isel(y=len(temp.y)//2,t=t).values.squeeze(),color="skyblue",alpha=0.5)        
        # plot depth line and scatter at centre y
        ax.plot(np.concatenate(([xmin],temp.x.values)),np.concatenate(([temp.isel(x=0,y=len(temp.y)//2,t=t).Watlev.item()*sc_wlv],-temp.isel(y=len(temp.y)//2).Botlev.values)),'k--')
        xr.plot.scatter(temp.isel(y=len(temp.y)//2,t=t,x=0),ax=ax,x="x",y="Watlev",c="k",s=200)
        if depmin: # runup line and scatter
            ax.plot([xrunup.isel(y=len(temp.y)//2,t=t).item()]*2,[zmin,temp.Ibp.isel(y=len(temp.y)//2,t=t).item()],c='r',ls='--')
            ax.scatter(xrunup.isel(y=len(temp.y)//2,t=t),temp.isel(y=len(temp.y)//2,t=t).Ibp,c="r",s=200)
            aax=ax.inset_axes(
                [0.7,0.2,0.2,0.2], transform=ax.transAxes)
            ibp=temp.Ibp.isel(y=len(temp.y)//2)
            if axis_off:
                [(i.set_xticklabels(""),i.set_yticklabels(""),i.axis([ibp.t.isel(t=0),ibp.t.isel(t=-1),ibpmin,ibpmax]),i.set_xlabel("t"),i.set_ylabel("z (SWL)"),i.grid("on"),i.axhline(y=0,c="k",ls="--",lw=2)) for i in [aax]]
            else:
                [(i.axis([ibp.t.isel(t=0),ibp.t.isel(t=-1),ibp.min().item(),ibp.max().item()]),i.set_xlabel("t [s]"),i.set_ylabel("z [m]"),i.grid("on")) for i in [aax]]
            aax.plot([ibp.isel(t=t).t.item()],[ibp.isel(t=t).item()],ls="",marker="o",c="r")
            aax.plot(ibp.isel(t=slice(0,t+1)).t.values,ibp.isel(t=slice(0,t+1)).values,ls="-",marker="",c="r")        
        # set axis properties
        if axis_off:
            ax.axis('off')
        else:
            # print time
            ax.text(0.5,1.1,'t = {:.1f} s'.format(temp.isel(t=t).t.item()),transform=ax.transAxes,ha='center',va='bottom')
            [(i.set_xlabel('X [m]'),i.set_ylabel('$\zeta,-d $ [m]')) for i in [ax]]
        [(i.set_xlim([xmin-(xmax-xmin)*.01,xmax]),i.set_ylim([zmin,max(zmax,vmax)])) for i in [ax]]
        ax.set_title("")
    def data_gen(ax,t=0): # lower 3D plot
        # create seabed colormap
        cmapb = 'YlOrBr_r'
        plt.cm.register_cmap(name='seabed',cmap=colors.LinearSegmentedColormap.from_list('seabed',
                        [(0,plt.cm.get_cmap(cmapb)(int(0.65*255))),
                        (1,plt.cm.get_cmap(cmapb)(int(.80*255)))]
                        ))
        # plot seabed
        (-temp.Botlev).plot.surface(ax=ax,rstride=1,cstride=1,cmap=plt.cm.get_cmap('seabed',(zmax-zmin)*4),vmax=zmax,vmin=zmin,alpha=0.9,zorder=zorder[0],add_colorbar=False) #seabed with increased arrays (to make it look nicer)        
        # select water level colormap
        colmap=plt.cm.get_cmap('Blues',256); colmap.set_bad(color='w',alpha=0) #set color for nan
        M=np.isnan(temp.Watlev.isel(t=t)) #set mask for nan and depths below hrun
        light = colors.LightSource(azim_light,elev_light) # azimuth and elevation; light source to be applied to water surface contour face color
        # plot water level
        (sc_wlv*temp.Watlev.isel(t=t)).plot.surface(ax=ax,rstride=1, cstride=1, linewidth=0,
                        antialiased=False, facecolors=light.shade(np.ma.masked_array(temp.Watlev.isel(t=t),mask=M),cmap=colmap,vmin=vmin*sc_wlv,vmax=vmax*sc_wlv),
                        zorder=zorder[1],alpha=.8,add_colorbar=False)
        # plot water level line and scatter at centre y and start of x, and profile of depth at centre y
        ax.plot(np.concatenate(([xmin],temp.x.values)),[temp.isel(y=len(temp.y)//2).y.values]*(len(temp.x)+1),np.concatenate(([temp.isel(x=0,y=len(temp.y)//2,t=t).Watlev.item()*sc_wlv],-temp.isel(y=len(temp.y)//2).Botlev.values)),'k--',zorder=zorder[2],alpha=1)
        ax.scatter(temp.isel(x=0).x,temp.isel(y=len(temp.y)//2).y,temp.isel(y=len(temp.y)//2,t=t,x=0).Watlev*sc_wlv,s=200,c='k')
        if depmin: # plot runup line (scatter not possible to see)
            ax.plot(xrunup.isel(t=t).values,xrunup.isel(t=t).y.values,temp.Ibp.isel(t=t).values*sc_wlv,'r-',zorder=zorder[3])
        # plot ymin and ymax planes
        for y in [0,-1]:
            ax.add_collection3d(Poly3DCollection([[(i.x.item(),i.y.item(),i.item()) for i in (-temp.Botlev).isel(y=y)]+\
                     [(i.x.item(),i.y.item(),zmin) for i in temp.isel(y=y).isel(x=slice(None,None,-len(temp.x)+1),t=0).Watlev]],
                                facecolor='burlywood'))
        # plot xmax plane
        for x in [-1]:
            ax.add_collection3d(Poly3DCollection([[(i.x.item(),i.y.item(),zmin) for i in temp.isel(x=x,y=slice(None,None,len(temp.y)-1),t=0).Watlev]+\
                     [(i.x.item(),i.y.item(),zmax) for i in temp.isel(x=x,y=slice(None,None,-len(temp.y)+1),t=0).Watlev]],
                                facecolor='burlywood'))            
        # set axis properties
        [(i.grid(False),i.set_title(""),i.set_box_aspect(aspect_ratio)) for i in [ax]]
        for i in [ax]:
            i.xaxis.pane.fill = False; i.yaxis.pane.fill = False; i.zaxis.pane.fill = False
            i.xaxis.labelpad = 10; i.yaxis.labelpad = 10
        if axis_off:
             ax.axis('off')
        else:
            [(i.set_xlabel('X [m]'),i.set_ylabel('Y [m]'),i.set_zlabel('Z [m]')) for i in [ax]]
        [(i.set_xlim([xmin-(xmax-xmin)*.01,xmax]),i.set_ylim([ymin,ymax]),i.set_zlim([zmin,max(zmax,vmax)])) for i in [ax]]
        ax.view_init(elev,azim) #elevation and azimuth
    
    def fig_rp(t=0):    
        fig = plt.figure(figsize=(15,10), constrained_layout=True)
        specc = fig.add_gridspec(ncols=1,nrows=2,height_ratios=[1.5,5])
        wlv_1d(fig.add_subplot(specc[0]),t=t)
        data_gen(fig.add_subplot(specc[1],projection='3d'),t=t)
        plt.savefig(f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png',dpi=dpi,bbox_inches='tight')
        plt.close()

    if not combine_only:
        if tpar==-1:
            print(f"Creating figures - total of {len(temp.t)} time steps")
            for t in range(len(temp.t)):
                if t%10==0: print(f"{t+1}/{len(temp.t)}")
                fig_rp(t=t)
            print(f"Creating gif")
            frames_gif(gf=len(temp.t),dt=dt,gif_name=gif_name,gif_path=gif_path,acc_factor=acc_factor)
            print(f"Deleting figures")
            delete_fig(gf=len(temp.t),gif_name=gif_name,gif_path=gif_path)
        else:
            #print(f"{tpar+1}/{len(temp.t)}") # for debugging
            number_pngs=len([i for i in os.listdir(f'{gif_path}{gif_name}/') if i.split('.')[-1]=='png'])
            if number_pngs%10==0: print(f"{number_pngs/len(temp.t)*100:.1f} %")
            fig_rp(t=tpar)       
    else:
            print(f"Creating gif")
            frames_gif(gf=len(temp.t),dt=dt,gif_name=gif_name,gif_path=gif_path,acc_factor=acc_factor)
            print(f"Deleting figures")
            delete_fig(gf=len(temp.t),gif_name=gif_name,gif_path=gif_path)


def frames_gif(gf,dt,gif_name,gif_path="",acc_factor=1):
    """Combine frames into gif file
    Args:
        gf (int): Number of frames.
        dt (float): time step (s).
        gif_name (str): Gif name
        gif_path (str, optional): path where to save gif. Default to "" (local path).
        acc_factor (float, optional): Ratio of frame duration by dt. Default to 1.
    """
    frames=[]
    for filename in [f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png'  for t in range(gf)]:
        frames.append(imageio.imread(filename))
    imageio.mimsave(f'{gif_path}{gif_name}.gif', frames, 'GIF', duration=dt/acc_factor)

def delete_fig(gf,gif_name,gif_path=""):
    """
    Delete original figures
    Args:
        gf (int): Number of frames.
        gif_name (str): Gif name
        gif_path (str, optional): path where to save gif. Default to "" (local path).
    """    
    for filename in [f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png'  for t in range(gf)]:
        os.remove(filename)
    os.rmdir(f'{gif_path}{gif_name}')


def create_2D_comp_gif(gif_name,ds,label,xmin=None,xmax=None,tmin=None,tmax=None,dt=None,zmin=None,zmax=None,depmin=None,acc_factor=1,gif_path="",axis_off=False,dpi=100):
    """
    Create 2D gif with water level animation for comparison between runs.
    Args:
        gif_name (str): Gif name
        ds (list): List with single xr data structure with 'x', 'Botlev' and 'Watlev'.
        label (list): List with names of run (for legend).
        xmin (float): minimum x-position (m). If None, ds.x.min().
        xmax (float): maximum x-position (m). If None, ds.x.max().
        tmin (float): minimum time (s). If None, ds.t.min().
        tmax (float): maximum time (s). If None, ds.t.max().
        dt (float): time step (s). If None, (ds.t.isel(t=1)-ds.t.isel(t=0)).
        zmin (float): minimum z-position (m). If None, (-ds.Botlev).min().
        zmax (float): maximum z-position (m). If None, ds.Watlev.max().
        acc_factor (float, optional): Ratio of frame duration by dt. Default to 1.
        gif_path (str, optional): path where to save gif. Default to "" (local path).
        axis_off (bool, optional): boolean for plotting axis. Default to False (i.e., axis is plotted).
        dpi (int): Image dpi. Default to 100.
    """
    warnings.filterwarnings('ignore')
    # Assign xmin, xmax, tmin, tmax, dt, zmin, and xmax if not defined
    xmin = xmin or min([i.x.min().item() for i in ds])
    xmax = xmax or min([i.x.max().item() for i in ds])
    tmin = tmin or min([i.t.min().item() for i in ds])
    tmax = tmax or min([i.t.max().item() for i in ds])
    dt = dt or (ds[0].t.isel(t=1)-ds[0].t.isel(t=0)).item()
    temp=[i.sel(x=slice(xmin,xmax),t=slice(tmin,tmax,math.ceil(dt/((i.t.isel(t=1)-i.t.isel(t=0)).values.item())))) for i in ds]
    zmin = zmin or min([(-i.Botlev).min().item() for i in temp])
    zmax = zmax or max([max(i.Watlev.max().item(),(-i.Botlev).max().item()) for i in temp])
    if not os.path.exists(gif_name): os.mkdir(f'{gif_path}{gif_name}')
    d=temp[0].Botlev # bottom level - only for first in the list
    wl=[i.Watlev for i in temp] # water level
    def fig_rp(t=0):
        _,ax=plt.subplots(constrained_layout=True,figsize=(10,3))
        ax=[ax]
        c=["k","r","b","g","gray"]
        # plot land and water level surfaces
        ax[0].fill_between(d.x,zmin+0*d.values.squeeze(),-d.values.squeeze(),color="peachpuff") # contour beach
        #ax[0].fill_between(d.x,-d.values.squeeze(),wl.isel(t=t).values.squeeze(),color="skyblue",alpha=0.5)
        [i.isel(t=t).plot(ax=ax[0],color=j,label=l) for i,j,l in zip(wl,c,label)]
        ax[0].legend(loc="lower right")
        ax[0].set_title("")
        # set axis properties
        if axis_off:
            ax[0].axis('off')
        else:
            [(i.set_xlabel('X [m]'),i.set_ylabel('$\zeta,-d $ [m]')) for i in ax]
            ax[0].plot([(-d).where(lambda x:x>=0,drop=True).isel(x=0).x.item(),xmax],[0,0],c="k",ls="--") # at SWL
            ax[0].text(0.5,1.1,'t = {:.1f} s'.format(temp[0].isel(t=t).t.item()),transform=ax[0].transAxes,ha='center',va='bottom') # print time      
        ax[0].axhline(c="k",ls="--") # at SWL
        ax[0].text(xmax,0,"SWL",ha="right")
        [(i.set_xlim([xmin,xmax]),i.set_ylim([zmin,zmax])) for i in ax]
        plt.savefig(f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png',dpi=dpi)
        plt.close() 

    print(f"Creating figures - total of {len(temp[0].t)} time steps")
    for t in range(len(temp[0].t)):
        if t%50==0: print(f"{t+1}/{len(temp[0].t)}")
        fig_rp(t=t)
        
    print(f"Creating gif")
    frames_gif(gf=len(temp[0].t),dt=dt,gif_name=gif_name,gif_path=gif_path)
    print(f"Deleting figures")
    delete_fig(gf=len(temp[0].t),gif_name=gif_name,gif_path=gif_path)

if __name__ == '__main__':
    pass

