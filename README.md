# SWASH post-processing

This project contains functions to post-process the results of the non-hyrostatic model [SWASH](https://swash.sourceforge.io/).

## To install this package, execute the following commands:
```
python -m pip install git+https://github.com/rfonsecadasilva/swash_pp
```
## Alternatively for developer mode:
```
git clone https://github.com/rfonsecadasilva/swash_pp.git
cd swash_pp
pip install -e .
```

### Example of usage for creating nc from SWASH mat:
```
from swash_pp import swash_mat2nc as snc
snc.mat2nc_mean_2D(path_run="~/run_swash/")
```

### Example of loading nc with mean output, instantaneous water levels and runup:
```
import xarray as xr
run ="examp_run":
ds=([xr.open_dataset(run+"mean_2D.nc"),xr.open_dataset(run+"ins_2D_WATL_INST.nc"),xr.open_dataset(run+"ins_2D_IBP_RUNUP.nc")])
ds[1]=ds[1].assign_coords(x=ds[0].x) # assuming mean and water levels at the same position, otherwise interpolate
ds=xr.merge([ds[0],ds[1],ds[-1].isel(x=0,y=0).drop_vars(['y','x'])])
ds
```

### Example of usage for creating gif animation with water level from xr data structure:
```
from swash_pp import swash_wl_gif as wlg
gif_name="examp"
zmin,zmax=-0.4,0.4
xmin,xmax=46,54
tmin,tmax,dt=60,70,3.25/20
wlg.create_gif(gif_name=gif_name,ds=ds,xmin=xmin,xmax=xmax,tmin=tmin,tmax=tmax,dt=dt,zmin=zmin,zmax=zmax)
```
