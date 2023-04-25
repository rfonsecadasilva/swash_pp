def xvector(ds):
    """Return list with variables that are x-direction vectors

    Args:
        ds (xr dataset): input xarray dataset

    Returns:
        list: list with vectors that are x-direction vectors
    """
    xvector=[]
    for varname in ds.data_vars:
        if varname[-1]=="x" or ds[varname].attrs and "x-direction" in ds[varname].attrs["long_name"] or "U-velocity" in ds[varname].attrs["long_name"]:
            xvector.append(varname)
    return xvector

def x_shoreline(ds):
    """Return x-coordinate of shoreline

    Args:
        ds (xr dataset): input xarray dataset with "Botlev" variable

    Returns:
        float: x-coordinate of shoreline offset
    """
    return ds.isel(y=0).where(lambda ds:ds["Botlev"]<=0,drop=True).x.isel(x=0).item()


def make_x_shoreline_off(ds,x_beach):
    """Transform x-coordinate of dataset from shoreline to offshore

    Args:
        ds (xr dataset): input xarray dataset
        x_beach (float): x-postion of beach (in m).

    Returns:
        ds (xr dataset): output xarray dataset
    """
    ds=ds.assign_coords(x=-(ds.x-x_beach))
    ds=ds.reindex(x=ds.x[::-1])
    for k in xvector(ds):
        ds[k]=-ds[k]
    return ds
    
if __name__ == '__main__':
    pass