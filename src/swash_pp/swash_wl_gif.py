import warnings, math, os
import matplotlib.pyplot as plt
import imageio

def create_gif(gif_name,ds,xmin,xmax,tmin,tmax,dt,zmin,zmax,gif_path=""):
    """
    Create gif with water level animation
    Args:
        gif_name (str): Gif name
        ds (xr data structure): Single data structure with Botlev, Watlev and Ibp.
        xmin (float): minimum x-position (m)
        xmax (float): maximum x-position (m)
        tmin (float): minimum time (s)
        tmax (float): maximum time (s)
        dt (float): time step (s)
        zmin (float): minimum z-position (m)
        zmax (float): maximum z-position (m)
        gif_path (str): path where to save gif. Default to "" (local path).
    """
    
    warnings.filterwarnings('ignore')
    if not os.path.exists(gif_name): os.mkdir(f'{gif_path}{gif_name}')
    temp=temp=ds.sel(x=slice(xmin,xmax),t=slice(tmin,tmax,math.ceil(dt/((ds.t.isel(t=1)-ds.t.isel(t=0)).values.item()))))
    d=temp.Botlev # bottom level
    wl=temp.Watlev # water level
    if "Ibp" in temp: ibp=temp.Ibp # instantanenous beach position

    def fig_rp(t=0):
        
        fig,ax=plt.subplots(constrained_layout=True,figsize=(10,3))
        ax=[ax]
        ax[0].fill_between(d.x,zmin+0*d.values.squeeze(),-d.values.squeeze(),color="peachpuff") # contour beach
        ax[0].axis('off')
        ax[0].fill_between(d.x,-d.values.squeeze(),wl.isel(t=t).values.squeeze(),color="skyblue")
        ax[0].axis([xmin,xmax,zmin,zmax])
        ax[0].axhline(c="k",ls="--")
        ax[0].text(xmax,0,"SWL",ha="right")
        if "Ibp" in temp:
            aax=ax[0].inset_axes(
                [0.7,0.2,0.2,0.2], transform=ax[0].transAxes)
            [(i.set_xticklabels(""),i.set_yticklabels(""),i.axis([ibp.t.isel(t=0),ibp.t.isel(t=-1),ibp.min().item(),ibp.max().item()]),i.set_xlabel("t"),i.set_ylabel("z (SWL)"),i.grid("on"),i.axhline(y=ibp.mean().item(),c="k",ls="--",lw=2)) for i in [aax]]
            aax.plot([ibp.isel(t=t).t.item()],[ibp.isel(t=t).item()],ls="",marker="o",c="r")
            aax.plot(ibp.isel(t=slice(0,t+1)).t.values,ibp.isel(t=slice(0,t+1)).values,ls="-",marker="",c="r")
        plt.savefig(f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png',dpi=100)
        plt.close() 

    def frames_gif(gf,dt):
        """
        Combine frames into gif file
        """
        
        frames=[]
        for filename in [f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png'  for t in range(gf)]:
            frames.append(imageio.imread(filename))
        imageio.mimsave(f'{gif_path}{gif_name}.gif', frames, 'GIF', duration=dt)          

    def delete_fig(gf):
        """
        Delete original figures
        """
        for filename in [f'{gif_path}{gif_name}/{gif_name}_Fig_{t:04d}.png'  for t in range(gf)]:
            os.remove(filename)
        os.rmdir(f'{gif_path}{gif_name}')
        
    print(f"Creating figures - total of {len(temp.t)} time steps")
    for t in range(len(temp.t)):
        if t%50==0: print(f"{t+1}/{len(temp.t)}")
        fig_rp(t=t)
        
    print(f"Creating gif")
    frames_gif(gf=len(temp.t),dt=dt)
    print(f"Deleting figures")
    delete_fig(gf=len(temp.t))


if __name__ == '__main__':
    pass