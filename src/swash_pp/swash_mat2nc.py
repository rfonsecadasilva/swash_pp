"""
Write xarray dataset from SWASH gridded mat and table output runs.
@author: rfonsecadasilva
"""

def swashdict(path_sc=""):
    """
    Return dictionaries with SWASH output metadata from source code.

    Args:
        path_sc (str, optional): path of SWASH source code if using developer mode (otherwise it parses SWASH 9.01 online version). Defaults to "".
    
    Returns:
        swash_dict (dict): dictionary with SWASH output metadata (keyword, name, and units; read from SwashInit.ftn90).
        attr_dict (dict): dictionary with SWASH version metadata.
        out_ind (dict): dictionary with ivtype (SWASH source code number for each keyword output) to type.
    """
    def swashunit(init):
        """
        Return dict swash_units with correspondance between ovunit and units (e.g., 'uh' corresponds to m).

        Args:
            init(list): list with SwashInit.ftn90 imported into python (with readlines()).

        Returns:
            swash_unit(dict): dict with correspondance between ovunit and units
        """
        swash_unit={}
        lr=[]
        lr.append([index for index,value in enumerate(init) if "! units" in value][0]+2)
        lr.append([index for index,value in enumerate(init) if "! some constants" in value][0]-1)
        for i in range(lr[0],lr[1]):
            swash_unit[init[i].split()[0]]=init[i].split("'")[-2]
        return swash_unit

    if path_sc!="": # read from local path
        init = open(path_sc+"SwashInit.ftn90","r").readlines()
    else: # read from SWASH official release
        import urllib.request
        import tarfile
        from io import BytesIO
        tmpfile = BytesIO()
        tmpfile.write(urllib.request.urlopen("https://swash.sourceforge.io/download/zip/swash-9.01.tar.gz").read())
        tmpfile.seek(0)
        tfile = tarfile.open(fileobj=tmpfile, mode="r:gz")
        init=tfile.extractfile('swash-9.01/SwashInit.ftn90').readlines()
        tmpfile.close()
        tfile.close()
        init=[i.decode('UTF-8') for i in init]
    
    swash_dict={}
    lr=0
    while [index for index,value in enumerate(init[lr:]) if "ivtype = " in value]:
        lr+= [index for index,value in enumerate(init[lr:]) if "ivtype = " in value][0]
        ent=[]
        ent.append(int(init[lr].strip().split("=")[1])) # ivtype (SWASH number for each output variable), e.g., 5 for BOTLEV
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovkeyw" in value][0] # Keyword, e.g., BOTL
        ent.append(init[lr].split("'")[-2])
        lr = [lr+ index for index,value in enumerate(init[lr:]) if "ovsnam" in value][0] # Short name, e.g., Botlev
        ent.append(init[lr].split("'")[-2])
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovlnam" in value][0] # Long name, e.g., Bottom level
        ent.append(init[lr].split("'")[-2])
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovunit" in value][0] # SWASH unit (to be converted), e.g., uh (for m)
        ent.append(init[lr].split()[-1])
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovsvty" in value][0] # Type of variable, e.g., 1 for scalar, 2 for angle, 3 vector, 4 tensor
        ent.append(int(init[lr].split()[-1]))        
        swash_dict[ent[2]]=[ent[1],ent[3],ent[4],ent[0],ent[5]] # key is short name, list has: 0 Keyword 1 Long name 2 SWASH Units 3 Ivtype 4 Type of var
    
    # translate SWASH untis (e.g., uh) to real units (e.g., m)
    swash_unit = swashunit(init)
    for value in swash_dict.values():
        if value[2] in swash_unit: value[2]=swash_unit[value[2]]

    attr_dict={"title": "SWASH output",
               "description": "SWASH output converted from series of mat files",
               "SWASH version": [i for i in init if "VERNUM" in i][0].strip().split("=")[1]
               }    

    # dict with ivtype (SWASH source code number for each keyword output) to type (sta or ins; no_grid, 2D, ke, or kc)
    # based on swash-9.01_further - https://github.com/rfonsecadasilva/swash-9.01_further
    # sta is static; ins is instantaneous; no_grid means constant output variables; 2D is twodimensional variable; ke is 3D defined at vertical cell edge, whereas kc is at cell center)
    out_ind={}
    out_ind["no_grid"]=[40,41] + list(range(101,116)) # no grid instantaneous
    out_ind["sta_2D"] = [1,2,3,5,21,22,23,24,25,26,33,34,35,36,37,42,43,44,46]+list(range(149,185))+[193,194]+list(range(201,216))+[217,218]
    out_ind["ins_2D"] = [4]+list(range(6,21))+list(range(27,33))+list(range(38,40))+[45]+list(range(195,200))+[216]+list(range(249,287))
    out_ind["sta_ke"] = [57,58,59,185,186]
    out_ind["ins_ke"] = [51,52,53,54,55,56,189,190]
    out_ind["sta_kc"] = [84,85,86,87,88,92,93,94,187,188]
    out_ind["ins_kc"] = list(range(71,84))+[89,90,91,191,192]

    # append dict for vector units
    for key,value in swash_dict.copy().items():
        if value[4]==3: #vector
            swash_dict["X-"+key]=["X-"+value[0]] + [value[1]+ " - x-component"] + value[2:]
            swash_dict["Y-"+key]=["Y-"+value[0]] + [value[1]+ " - y-component"] + value[2:]
    return swash_dict,attr_dict,out_ind

def load_grid_req(path_run,run_file="run.sws"):
    """
    Parse SWASH input file and return relevant information for creating nc files.

    Args:
        path_run (str): path where run file (*.sws) is located.
        run_file (str, optional): run file name. Defaults to "run.sws".
    
    Returns:
        frame (list): list with grids (requested with BLOCK command) - id, y, and x
        reqn (list): list with gridded requests - id, type (mean or inst), and list of requested variables
    """
    import numpy as np
    sws=open(path_run+run_file,"r").readlines()
    # Make list with requests (id, type, fname and list of variables)
    req= [i for i in sws if i[:3]=='BLO' and '.mat' in i]
    # req: name, type (mean or inst), list of requested variables
    reqn = [ [i.split()[1],
            ['mean','inst']['OUT' in i], # mean and instantaneous variables
            [j.replace("'","") for j in i.split() if ".mat" in j][0], # name of .mat file
            [j for j in [i.split("OUT")[0].strip().split() for i in [i.split(".mat' ")[1] ]][0] if not j.isnumeric and "LAY" not in j]] # list of requested variables
            for i in req]
    # Make list with grids (id, y and x)
    fname=[i.split()[1] for i in sws if i[:3]=='FRA'] #frame names
    frame={} # id, y and x
    for f in fname:
        frame[f]=[[np.linspace(float(i.split()[3]),float(i.split()[3])+float(i.split()[6]),int(i.split()[8])+1) for i  in sws if i[:3]=='FRA' and  f==i.split()[1]][0], #y
                            [np.linspace(float(i.split()[2]),float(i.split()[2])+float(i.split()[5]),int(i.split()[7])+1) for i  in sws if i[:3]=='FRA' and  f==i.split()[1]][0]] #x
    # extract compgrid if regular grid
    for i in sws:
        if i[:5]=='CGRID' and "CURV" not in i and "UNSTRUC" not in i:
            first_index=1
            if "REG" in i: first_index+=1
            cgrid={"'COMPGRID'":[np.linspace(float(i.split()[first_index+1]),float(i.split()[first_index+1])+float(i.split()[first_index+4]),int(i.split()[first_index+6])+1), #y
                                 np.linspace(float(i.split()[first_index]),float(i.split()[first_index])+float(i.split()[first_index+3]),int(i.split()[first_index+5])+1)]
                  }
            frame={**frame,**cgrid}
            break
    return frame,reqn

def load_table_req(path_run,run_file="run.sws"):
    """
    Parse SWASH input file and return relevant information for creating nc files.

    Args:
        path_run (str): path where run file (*.sws) is located.
        run_file (str, optional): run file name. Defaults to "run.sws".
    
    Returns:
        reqtable (list): list with table  requests - id, type (mean or inst), and list of requested variables.
        point (list): list with table points (requested with POI command) - id, y, and x.
        nV (int): number of vertical layers.
    """
    import numpy as np
    sws=[i for i in open(path_run+run_file,"r").readlines() if i[0]!="$"]
    nV=max(list(set([int(i.split()[1]) if i[:3]=='VER' else 1 for i in sws]))) # number of vertical layers
    mode = sorted(list(set(["ONED" if i[:]=='MODE' and "ONED" in i else "TWOD" for i in sws])))[0] # extract mode
    # Make list with table requests (id, type, fname and list of variables)
    reqt= [i for i in sws if i[:3]=='TAB' in i]         
     # reqt: name, type (mean or inst), list of requested variables
    reqtable = [ [i.split()[1],
            ['mean','inst']['OUT' in i], # mean and instantaneous variables
            [j.replace("'","") for j in i.split() if i.split(".")[1].split()[0][:-1] in j][0], # name of table file
            [k for k in i.split("OUT")[0].split(".")[1].split()[1:]]] # list of requested variables 
            for i in reqt]
    # make list with table
    fname=[i.split()[1] for i in sws if i[:3]=='POI' or i[:3]=='FRA'] #point names (note that table can also use grids)
    point={} # id, y and x
    for f in fname:
        for i in sws:
            if (i[:3]=='POI') and (f==i.split()[1]):
                if i.split()[2][0].isdigit():
                    point[f]=[[float(j) for j in i.split()[3:][::2]],[float(j) for j in i.split()[2:][::2]]]
                else:
                    with open(path_run+i.split()[3][1:-1],"r") as file:
                        tab=file.readlines()
                        point[f]=[[float(j.strip().split()[1]) for j in tab],[float(j.strip().split()[0]) for j in tab]]
            elif (i[:3]=='FRA') and (f==i.split()[1]):
                point[f]=[np.linspace(float(i.split()[3]),float(i.split()[3])+float(i.split()[6]),int(i.split()[8])+1), #y,
                          np.linspace(float(i.split()[2]),float(i.split()[2])+float(i.split()[5]),int(i.split()[7])+1)]
    # extract compgrid if regular grid
    for i in sws:
        if i[:5]=='CGRID' and "CURV" not in i and "UNSTRUC" not in i:
            first_index=1
            if "REG" in i: first_index+=1
            cgrid={"'COMPGRID'":[np.linspace(float(i.split()[first_index+1]),float(i.split()[first_index+1])+float(i.split()[first_index+4]),int(i.split()[first_index+6])+1), #y
                                 np.linspace(float(i.split()[first_index]),float(i.split()[first_index])+float(i.split()[first_index+3]),int(i.split()[first_index+5])+1)]
                  }
            point={**point,**cgrid}
            break
    return reqtable,point,nV

def mat2nc(path_run,path_sc="",run_file="run.sws",save_nc=False):
    """
    Convert all gridded matlab output into xr dataset and return list with all of them.

    Args:
        path_run (str): path where run file (*.sws) is located.
        path_sc (str, optional): path of SWASH source code if using developer mode (otherwise it parses SWASH 9.01 online version). Defaults to "".
        run_file (str, optional): run file name. Defaults to "run.sws".
        save_nc (bool, optional): option to save each nc file. Defaults to True.
    
    Returns:
        ds_out (list): list with xarray datasets for each SWASH matlab output.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    import re,os
    frame,reqn=load_grid_req(path_run,run_file=run_file)
    swash_dict,attr_dict,out_ind=swashdict(path_sc=path_sc)
    ds_out=[]
    for req in reqn:
        req[-1] = list(map(lambda x: x.replace('HSIG', 'HS'), req[-1]))
        if os.path.exists(path_run+req[2]): # check whether matlab file exists
            if all([i in [j[0] for j in swash_dict.values()] for i in req[-1]]): # check whether all output keywords can be identified
                out={}
                out.update(sio.loadmat(path_run+req[2])) # load swash output mat file
                # for vector quantities, create X- and Y- at the start of keys (consistent with swash_dict)
                out={(i if all([j not in i for j in ["_x","_y"]]) else "X-"+i.replace("_x","") if "_x" in i else "Y-"+i.replace("_y","")):j  for i,j in out.items()}
                # filter based on static/instantaneous 2D/3D quantities - iterate over output request for each grid
                out_sta_2D,out_ins_2D,out_sta_3D,out_ins_3D=\
                    [dict(filter(lambda item:swash_dict.get(re.split('(\d+)',item[0])[0].strip("_"),[""]*4)[3] in i,out.items()))
                            for i in [out_ind["sta_2D"],out_ind["ins_2D"],out_ind["sta_ke"]+out_ind["sta_kc"],out_ind["ins_ke"]+out_ind["ins_kc"]]]
                ds_sta_2D,ds_ins_2D,ds_sta_3D,ds_ins_3D="","","",""
                if out_sta_2D:
                    ds_sta_2D=xr.Dataset(
                                data_vars={key: (("y", "x"), value[::-1,:],{"standard_name": key,
                                                                            "long_name": swash_dict[key][1],
                                                                            "units":swash_dict[key][2]})
                                        for key,value in out_sta_2D.items() if key not in ('Xp','Yp')},
                                        attrs=attr_dict
                                )
                if out_ins_2D:
                    ds_ins_2D=[]
                    for var in list(set([re.split('(\d+)',i)[0].strip("_") for i in out_ins_2D.keys()])): #name of variable 
                        temp=dict(filter(lambda item:item[0].split("_")[0] == var,out_ins_2D.items()))
                        ds_ins_2D.append(xr.Dataset(
                                    data_vars={var: (("y", "x","t"), (np.stack(list(temp.values()), axis=-1))[::-1],
                                                                        {"standard_name": var,
                                                                        "long_name": swash_dict[var][1],
                                                                        "units": swash_dict[var][2]}),
                                            },
                                    coords={"t": [float(key[-10:-8])*60*60 +float(key[-8:-6])*60 +float(key[-6:-4]) + float(key[-3:])/1000 for key in temp.keys()]
                                        },
                                        attrs=attr_dict
                                    ))
                    ds_ins_2D=xr.merge(ds_ins_2D)
                if out_sta_3D:
                    ds_sta_3D=[]
                    for var in list(set([re.split('(\d+)',i)[0].strip("_") for i in out_sta_3D.keys()])): #name of variable
                        temp=dict(filter(lambda item:re.split('(\d+)',item[0])[0] == var,out_sta_3D.items()))
                        okv_k=[int(re.split('(\d+)',i)[1]) for i in temp.keys()]# k axis
                        ds_sta_3D.append(xr.Dataset(
                                    data_vars={var: (("y", "x",["kc","ke"][okv_k[0]==0]), np.stack(list(temp.values()),axis=-1)[::-1],
                                                                                {"standard_name":var,
                                                                                "long_name": swash_dict[var][1],
                                                                                "units": swash_dict[var][2]})},
                                    coords={["kc","ke"][okv_k[0]==0]: okv_k},
                                            attrs=attr_dict
                                        ))
                    ds_sta_3D=xr.merge(ds_sta_3D)
                if out_ins_3D:
                    ds_ins_3D=[]
                    for var in list(set([re.split('(\d+)',i)[0].strip("_") for i in out_ins_3D.keys()])): #name of variable
                        temp=dict(filter(lambda item:re.split('(\d+)',item[0])[0] == var,out_ins_3D.items()))
                        okv_k=list(set([int(re.split('(\d+)',i)[1]) for i in temp.keys()])) # k axis
                        temp=[dict(filter(lambda i:int(i[0].split("_")[1][1:])==k,temp.items())) for k in okv_k]# filter by layer k
                        temp_t=temp.copy() # for time
                        temp=[np.stack(list(i.values()), axis=-1)[::-1] for i in temp] # time stack
                        temp= np.stack([v for v in temp],axis=-2) #layer stack
                        ds_ins_3D.append(xr.Dataset(
                                    data_vars={var: (("y", "x",["kc","ke"][okv_k[0]==0],"t"), temp,
                                                                                {"standard_name":var,
                                                                                "long_name": swash_dict[var][1],
                                                                                "units": swash_dict[var][2]})},
                                    coords={["kc","ke"][okv_k[0]==0]: okv_k,
                                            "t":[float(key[-10:-8])*60*60 +float(key[-8:-6])*60 +float(key[-6:-4]) + float(key[-3:])/1000 for key in temp_t[0].keys()]},
                                            attrs=attr_dict
                                        ))
                    ds_ins_3D=xr.merge(ds_ins_3D)             
                ds_m=[j for i,j in zip([out_sta_2D,out_ins_2D,out_sta_3D,out_ins_3D],[ds_sta_2D,ds_ins_2D,ds_sta_3D,ds_ins_3D]) if i]
                ds=xr.merge(ds_m)
                if 'kc' in ds:
                    ds.kc.attrs = {"standard_name": 'kc',"long_name": 'Vertical axis (cell centre)', "units": '',"axis":"kc"}
                if 'ke' in ds:
                    ds.ke.attrs = {"standard_name": 'ke',"long_name": 'Vertical axis (cell edge)', "units": '',"axis":"ke"}
                if 't' in ds:
                    ds.t.attrs = {"standard_name": 'Tsec',"long_name": swash_dict['Tsec'][1], "units": swash_dict['Tsec'][2],"axis":"t"}        
                if req[0] in frame.keys(): # for regular grids
                    ds=ds.assign_coords(y=frame[req[0]][0])
                    ds.y.attrs = {"standard_name": "Yp","long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
                    ds=ds.assign_coords(x=frame[req[0]][1])
                    ds.x.attrs = {"standard_name": "Xp","long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}         
                else:
                    if 'Yp' in out:
                        ds=ds.assign_coords(y=out['Yp'][:,0][::-1])
                        ds.y.attrs = {"standard_name": "Yp","long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
                    if 'Xp' in out:
                        ds=ds.assign_coords(x=out['Xp'][0,:])
                        ds.x.attrs = {"standard_name": "Xp","long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}         
                if save_nc:
                    print("Saving "+path_run+f"{req[2][:-4]}.nc")
                    ds.to_netcdf(path_run+f"{req[2][:-4]}.nc")
                ds_out.append(ds)
            else:
                print(f"Warning: output keyword {[i for i in req[-1] if i not in [j[0] for j in swash_dict.values()]]} does not valid metadata and {path_run+req[2]} has been skipped")
        else:
            print(f"Warning: {path_run+req[2]} does not exist and has been skipped")
            continue

    return ds_out

def tab2nc(path_run,path_sc="",run_file="run.sws",save_nc=False):
    """
    Convert all table output into xr dataset and return list with all of them.

    Args:
        path_run (str): path where run file (*.sws) is located.
        path_sc (str, optional): path of SWASH source code if using developer mode (otherwise it parses SWASH 9.01 online version). Defaults to "".
        run_file (str, optional): run file name. Defaults to "run.sws".
        save_nc (bool, optional): option to save each nc file. Defaults to True.
    
    Returns:
        ds_out (list): list with xarray datasets for each SWASH table output.    
    """
    import xarray as xr
    import numpy as np
    import os
    def custom_index(first_index,chunk,skip,length):
        """Retrieve custom indexing.

        Args:
            first_index (int): first index of np array.
            chunk (int): chunk length.
            skip (int): skip.
            length (int): length of np array.

        Returns:
            out (np array): np array with custom indexes.
        """
        return np.concatenate([range(j,j+chunk) for j in range(first_index,length,skip)])

    reqtable,point,nV=load_table_req(path_run,run_file=run_file)
    swash_dict,attr_dict,out_ind=swashdict(path_sc=path_sc)
    swash_dict={value[0]:[key]+value[1:] for key,value in swash_dict.items()}
    ds_out=[]
    for tab in reqtable: #name, of grid type (mean or inst), name of tbl, list of requested variables
        if os.path.exists(path_run+tab[2]): # check whether table file exists
            if all([i in swash_dict.keys() for i in tab[-1]]): # check whether all output keywords can be identified
                # modify tab to account for vectors
                tab_rev=[]
                [tab_rev.append(i) if swash_dict[i][4] <3 else (tab_rev.append("X-"+i),tab_rev.append("Y-"+i))  for i in tab[-1]]
                # initiate list with xr datasets
                ds=[]
                # load file
                file=[i.strip() for i in open(path_run+tab[2]).readlines()]
                if file[0][0]=="S": # condition for files starting with "SWASH   1" ..." (most likely files with at least one variable with vertical components)
                    file=file[6+4*(len(tab_rev)):] 
                else: # condition for ignoring lines starting with %
                    temp=0
                    for i in range(len(file)):
                        if file[i][0]!="%":
                            temp=i; break
                    file=file[temp:]
                file=np.concatenate([[float(k) for k in j.split()] for j in file]) # one line with all values
                # 2D files are output at same line whereas if any 3D each output is given at a different line
                ndim=[nV+1 if swash_dict[i][3] in out_ind["sta_ke"]+out_ind["ins_ke"] else nV if swash_dict[i][3] in out_ind["sta_kc"]+out_ind["ins_kc"] else 1 for i in tab_rev]
                # initiate dict that will contain the whole dataset
                out={}
                for i in range(len(tab_rev)): # iterate over each keyword requested
                    out[swash_dict[tab_rev[i]][0]]=file[custom_index(sum(ndim[:i]),ndim[i],int(np.dot(ndim,np.ones(len(ndim)))),len(file))] # retrieve values for keyword
                    if ndim[i]==1: # reshape into station and time
                        if tab[0] != "'NOGRID'":
                            out[swash_dict[tab_rev[i]][0]]=np.vstack([out[swash_dict[tab_rev[i]][0]][j::len(point[tab[0]][1])] for j in range(len(point[tab[0]][1]))])  # reorder to groups of stations (rather than default groups of time) and reshape into station and time
                            if swash_dict[tab_rev[i]][0] not in ["Time","Tsec"]: # exclude time as it consists of a coordinate (see below)
                                ds.append(xr.Dataset(
                                    data_vars={swash_dict[tab_rev[i]][0]: (("stations", "t"), out[swash_dict[tab_rev[i]][0]],
                                                                                {"standard_name": swash_dict[tab_rev[i]][0],
                                                                                "long_name": swash_dict[tab_rev[i]][1],
                                                                                "units": swash_dict[tab_rev[i]][2]})
                                            },
                                            attrs=attr_dict
                                            ))
                        else:
                            if swash_dict[tab_rev[i]][0] not in ["Time","Tsec"]:
                                ds.append(xr.Dataset(
                                    data_vars={swash_dict[tab_rev[i]][0]: (("t"), out[swash_dict[tab_rev[i]][0]],
                                                                                {"standard_name": swash_dict[tab_rev[i]][0],
                                                                                "long_name": swash_dict[tab_rev[i]][1],
                                                                                "units": swash_dict[tab_rev[i]][2]})
                                            },
                                            attrs=attr_dict
                                            ))                            
                    else: # reshape into station, time and k
                        out[swash_dict[tab_rev[i]][0]]=np.array([[out[swash_dict[tab_rev[i]][0]][j*ndim[i]+k::len(point[tab[0]][1])*ndim[i]] for j in range(len(point[tab[0]][1]))] for k in range(ndim[i])]) # reorder and reshape to k,station,time
                        out[swash_dict[tab_rev[i]][0]]=np.transpose(out[swash_dict[tab_rev[i]][0]],(1,2,0)) # transpose to station,time,k
                        ds.append(xr.Dataset(
                                    data_vars={swash_dict[tab_rev[i]][0]: (("stations","t",["kc","ke"][swash_dict[tab_rev[i]][3] in out_ind["sta_ke"]+out_ind["ins_ke"]]), out[swash_dict[tab_rev[i]][0]],
                                                                                {"standard_name": swash_dict[tab_rev[i]][0],
                                                                                "long_name": swash_dict[tab_rev[i]][1],
                                                                                "units": swash_dict[tab_rev[i]][2]})
                                            },
                                            attrs=attr_dict
                                        ))
                ds=xr.merge(ds)
                if 'Time' in out.keys():
                    t="Time"
                    if tab[0] != "'NOGRID'": ds=ds.assign_coords(t=out[t][0,:])
                    else: ds=ds.assign_coords(t=out[t][:])
                    ds.t.attrs = {"standard_name": t,"long_name": swash_dict['TIME'][1], "units": swash_dict["TIME"][2],"axis":"t"}
                elif 'Tsec' in out:
                    t="Tsec"
                    if tab[0] != "'NOGRID'": ds=ds.assign_coords(t=out[t][0,:])
                    else: ds=ds.assign_coords(t=out[t][:])
                    ds.t.attrs = {"standard_name": t,"long_name": swash_dict["TSEC"][1], "units": swash_dict["TSEC"][2],"axis":"t"}
                if tab[0] != "'NOGRID'":
                    if len(point[tab[0]][0])==len(point[tab[0]][1]):
                        ds['station_y']=(("stations"),point[tab[0]][0])
                    else:
                        ds['station_y']=(("stations"),np.repeat(point[tab[0]][0][0],len(point[tab[0]][1])))
                    ds['station_x']=(("stations"),point[tab[0]][1])
                    ds.station_y.attrs = {"standard_name": swash_dict['YP'][0],"long_name": swash_dict['YP'][1], "units": swash_dict['YP'][2],"axis":"Y"}
                    ds.station_x.attrs = {"standard_name": swash_dict['XP'][0],"long_name": swash_dict['XP'][1], "units": swash_dict['XP'][2],"axis":"X"}    
                if save_nc:
                    print("Saving "+path_run+f"{tab[2][:-4]}.nc")
                    ds.to_netcdf(path_run+f"{tab[2][:-4]}.nc")
                ds_out.append(ds)
            else:
                print(f"Warning: output keyword {[i for i in tab[-1] if i not in swash_dict.keys()]} does not valid metadata and {path_run+tab[2]} has been skipped")
        else:
            print(f"Warning: {path_run+tab[2]} does not exist and has been skipped")
            continue
    return ds_out

if __name__ == '__main__':
    pass