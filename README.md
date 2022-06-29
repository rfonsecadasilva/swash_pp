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

### Example of usage:
```
from swash_pp import swash_mat2nc as snc
snc.mat2nc_mean_2D(path_run="~/run_swash/")
```