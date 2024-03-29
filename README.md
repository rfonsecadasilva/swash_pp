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

### Example of usage for creating list of xarray datasets from SWASH mat files:
```
from swash_pp import swash_mat2nc as snc
ds=snc.mat2nc(path_run="~/run_swash/",run_file="run.sws")
```

### Example of usage for creating gif animation with water level from xarray data structure ds including 'x', 'Botlev', and 'Watlev':
```
from swash_pp import swash_wl_gif as wlg
gif_name="examp"
wlg.create_2D_gif(gif_name=gif_name,ds=ds)
```
