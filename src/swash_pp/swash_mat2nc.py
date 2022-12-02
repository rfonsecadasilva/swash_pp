"""
Created on May 04 2021
Write nc files from SWASH gridded mat output runs.
@author: rfonsecadasilva
"""

def swashdict(path_sc=""):
    """
    Return swash_dict (dictionary with output metadata)
    path_sc is the path of SWASH source code if using developer mode.
    """
    def swashunit(init):
        """
        Return dict swash_units (e.g., uh corresponds to m)
        """
        swash_unit={}
        lr=[]
        lr.append([index for index,value in enumerate(init) if "! units" in value][0]+2)
        lr.append([index for index,value in enumerate(init) if "! some constants" in value][0]-1)
        for i in range(lr[0],lr[1]):
            swash_unit[init[i].split()[0]]=init[i].split("'")[-2]
        return swash_unit

    if path_sc!="":
        init = open(path_sc+"SwashInit.ftn90","r").readlines()
    else:
        import urllib.request
        import tarfile
        from io import BytesIO
        tmpfile = BytesIO()
        tmpfile.write(urllib.request.urlopen("https://swash.sourceforge.io/download/zip/swash-7.01.tar.gz").read())
        tmpfile.seek(0)
        tfile = tarfile.open(fileobj=tmpfile, mode="r:gz")
        init=tfile.extractfile('swash/SwashInit.ftn90').readlines()
        tmpfile.close()
        tfile.close()
        init=[i.decode('UTF-8') for i in init]
    
    swash_dict={}
    lr=0
    while [index for index,value in enumerate(init[lr:]) if "ivtype = " in value]:
        lr+= [index for index,value in enumerate(init[lr:]) if "ivtype = " in value][0]
        ent=[]
        ent.append(int(init[lr].strip().split("=")[1])) # ivtype
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovkeyw" in value][0]
        ent.append(init[lr].split("'")[-2])
        lr = [lr+ index for index,value in enumerate(init[lr:]) if "ovsnam" in value][0]
        ent.append(init[lr].split("'")[-2])
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovlnam" in value][0]
        ent.append(init[lr].split("'")[-2])
        lr = [lr+index for index,value in enumerate(init[lr:]) if "ovunit" in value][0]
        ent.append(init[lr].split()[-1])
        swash_dict[ent[1]]=[ent[1],ent[3],ent[4],ent[2],ent[0]]
    
    swash_unit = swashunit(init)
    for value in swash_dict.values():
        if value[2] in swash_unit: value[2]=swash_unit[value[2]] # better than get, as string sometimes does not translate into acronym

    attr_dict={"title": "SWASH output",
               "description": "SWASH output converted from series of mat files",
               "SWASH version": [i for i in init if "VERNUM" in i][0].strip().split("=")[1]
               }    

    return swash_dict,attr_dict

def load_req(path_run,run_file="run.sws"):
    """
    Return frame (list with grids - id, y and x) and reqn (list with request info - name, type (mean or inst), list of requested variables)
    """
    import numpy as np
    sws=open(path_run+run_file,"r").readlines()
    # check nV
    nV=int([i.split()[1] for i in sws if i[:3]=='VER'][0])
    # Make list with requests (id, type, fname and list of variables)
    req= [i for i in sws if i[:3]=='BLO' and '.mat' in i]
    # req: name, type (mean or inst), list of requested variables
    reqn = [ [i.split()[1],
            ['mean','inst']['OUT' in i], # mean and instantaneous variables
            [j.replace("'","") for j in i.split() if ".mat" in j][0], # name of .mat file
            [i.split("OUT")[0].strip().split() for i in [i.split(".mat' ")[1] ]][0]] # list of requested variables
            for i in req]
    # Make list with table requests (id, type, fname and list of variables)
    reqt= [i for i in sws if i[:3]=='TAB' in i]            
     # reqt: name, type (mean or inst), list of requested variables
    reqtable = [ [i.split()[1],
            'inst', # assuming instantaneous variables
            [j.replace("'","") for j in i.split() if ".tbl" in j][0], # name of .tbl file
            [i.split("OUT")[0].strip().split() for i in [i.split(".tbl' ")[1] ]][0]] # list of requested variables
            for i in reqt]
    # Make list with grids (id, y and x)
    fname=[i.split()[1] for i in sws if i[:3]=='FRA'] #frame names
    frame={} # id, y and x
    for f in fname:
        frame[f]=[[np.linspace(float(i.split()[3]),float(i.split()[3])+float(i.split()[6]),int(i.split()[8])+1) for i  in sws if i[:3]=='FRA' and  f==i.split()[1]][0], #y
                            [np.linspace(float(i.split()[2]),float(i.split()[2])+float(i.split()[5]),int(i.split()[7])+1) for i  in sws if i[:3]=='FRA' and  f==i.split()[1]][0]] #x
    # make list with table
    fname=[i.split()[1] for i in sws if i[:3]=='POI'] #point names
    point={} # id, y and x
    for f in fname:
        for i in sws:
            if i[:3]=='POI' and f==i.split()[1]:
                if i.split()[2][0].isdigit():
                    point[f]=[[float(j) for j in i.split()[3:][::2]],[float(j) for j in i.split()[2:][::2]]]
                else:
                    with open(path_run+i.split()[2],"r").readlines() as tab:
                        point[f]=[[j.strip().split()[1] for j in tab],[j.strip().split()[0] for j in tab]]
    # cgrid for regular grid
    cgrid={"'COMPGRID'":[np.linspace(float(i.split()[3]),float(i.split()[3])+float(i.split()[6]),int(i.split()[8])+1), #y
                            np.linspace(float(i.split()[2]),float(i.split()[2])+float(i.split()[5]),int(i.split()[7])+1)]
                for i in sws if i[:5]=='CGRID' and i.split()[1][:3]=='REG'}
    # cgrid for regular grid
    cgrid={"'COMPGRID'":[np.linspace(float(i.split()[3]),float(i.split()[3])+float(i.split()[6]),int(i.split()[8])+1), #y
                            np.linspace(float(i.split()[2]),float(i.split()[2])+float(i.split()[5]),int(i.split()[7])+1)]
                for i in sws if i[:5]=='CGRID' and i.split()[1][:3]=='REG'}
    if cgrid: frame={**frame,**cgrid}
    else:
        import scipy.io as sio
        for i in reqn:
            if 'YP' and 'XP' in i[3]:
                out={}
                out.update(sio.loadmat(path_run+i[2]))
                out= dict(filter(lambda item:item[0] in ('Xp','Yp'),out.items())) # basic filter
                cgrid={"'COMPGRID'":[list(out['Yp'][:,0][::-1]),list(out['Xp'][0,:])]}
                frame={**frame,**cgrid}
    return frame,reqn,nV,reqtable,point

def mat2nc_all(path_run,path_sc="",run_file="run.sws"):
    """
    Save one nc for each gridded matlab output
    path_sc is the path of SWASH source code if using developer mode.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    import re
    frame,reqn,nV,reqtable,point=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    # ivtype to type - based on swash-8.01_further
    #no_grid=list(range(101,114)) # no grid instantaneous
    sta_2D = [1,2,3,5,21,22,23,24,25,26,33,34,35,36,37,42,43,44,46]+list(range(149,185))+[193,194]+list(range(201,216))+[217,218]
    ins_2D = [4]+list(range(6,21))+list(range(27,33))+list(range(38,42))+[45]+list(range(195,200))+[216]+list(range(249,287))
    sta_ke = [57,58,59,185,186]
    ins_ke = [51,52,53,54,55,56,189,190]
    sta_kc = [84,85,86,87,88,92,93,94,187,188]
    ins_kc = list(range(71,84))+[89,90,91,191,192]
    for req in reqn:
        out={}
        out.update(sio.loadmat(path_run+req[2])) # load mat file
        req[-1] = list(map(lambda x: x.replace('HSIG', 'HS'), req[-1]))
        # filter based on static/instantaneous 2D/3D quantities - iterate over output request for each grid
        out_sta_2D,out_ins_2D,out_sta_3D,out_ins_3D=\
            [dict(filter(lambda item:re.split('(\d+)',item[0])[0].strip("_") in ([swash_dict[var][3] for var in req[-1] if swash_dict[var][4] in i]),out.items()))
                    for i in [sta_2D,ins_2D,sta_ke+sta_kc,ins_ke+ins_kc]]
        ds_sta_2D,ds_ins_2D,ds_sta_3D,ds_ins_3D="","","",""
        if out_sta_2D:
            ds_sta_2D=xr.Dataset(
                        data_vars={key: (("y", "x"), value[::-1,:],{"standard_name": key,
                                                                    "long_name": swash_dict.get(key,['']*2)[1],
                                                                    "units": swash_dict.get(key,['']*3)[2]})
                                for key,value in out_sta_2D.items() if key not in ('Xp','Yp')},
                                attrs=attr_dict
                        )
        if out_ins_2D:
            ds_ins_2D=[]
            for var in [i for i in req[-1] if swash_dict[i][4] in ins_2D]: #name of variable 
                temp=dict(filter(lambda item:item[0].split("_")[0] in swash_dict[var][3],out_ins_2D.items()))
                ds_ins_2D.append(xr.Dataset(
                            data_vars={swash_dict[var][3]: (("y", "x","t"), (np.stack(list(temp.values()), axis=-1))[::-1],
                                                                {"standard_name": swash_dict[var][3],
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
            for var in [i for i in req[-1] if swash_dict[i][4] in sta_ke+sta_kc]: #name of variable
                temp=dict(filter(lambda item:re.split('(\d+)',item[0])[0] in swash_dict[var][3],out_sta_3D.items()))
                okv_k=[int(re.split('(\d+)',i)[1]) for i in temp.keys()]# k axis
                ds_sta_3D.append(xr.Dataset(
                            data_vars={swash_dict[var][3]: (("y", "x",["kc","ke"][okv_k[0]==0]), np.stack(temp.values(),axis=-1)[::-1],
                                                                        {"standard_name":swash_dict[var][3],
                                                                        "long_name": swash_dict[var][1],
                                                                        "units": swash_dict[var][2]})},
                            coords={["kc","ke"][okv_k[0]==0]: okv_k},
                                    attrs=attr_dict
                                ))
            ds_sta_3D=xr.merge(ds_sta_3D)
        if out_ins_3D:
            ds_ins_3D=[]
            for var in [i for i in req[-1] if swash_dict[i][4] in ins_ke+ins_kc]: #name of variable
                temp=dict(filter(lambda item:re.split('(\d+)',item[0])[0] in swash_dict[var][3],out_ins_3D.items()))
                okv_k=list(set([int(re.split('(\d+)',i)[1]) for i in temp.keys()])) # k axis
                temp=[dict(filter(lambda i:int(i[0].split("_")[1][1:])==k,temp.items())) for k in okv_k]# filter by layer k
                temp_t=temp.copy() # for time
                temp=[np.stack(list(i.values()), axis=-1)[::-1] for i in temp] # time stack
                temp= np.stack([v for v in temp],axis=-2) #layer stack
                ds_ins_3D.append(xr.Dataset(
                            data_vars={swash_dict[var][3]: (("y", "x",["kc","ke"][okv_k[0]==0],"t"), temp,
                                                                        {"standard_name":swash_dict[var][3],
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
            ds.t.attrs = {"standard_name": swash_dict['TSEC'][0],"long_name": swash_dict['TSEC'][1], "units": swash_dict['TSEC'][2],"axis":"t"}        
        if req[0] in frame.keys(): # for regular grids
            ds=ds.assign_coords(y=frame[req[0]][0])
            ds.y.attrs = {"standard_name": swash_dict['YP'][0],"long_name": swash_dict['YP'][1], "units": swash_dict['YP'][2],"axis":"Y"}
            ds=ds.assign_coords(x=frame[req[0]][1])
            ds.x.attrs = {"standard_name": swash_dict['XP'][0],"long_name": swash_dict['XP'][1], "units": swash_dict['XP'][2],"axis":"X"}         
        else:
            if 'Yp' in out:
                ds=ds.assign_coords(y=out['Yp'][:,0][::-1])
                ds.y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
            if 'Xp' in out:
                ds=ds.assign_coords(x=out['Xp'][0,:])
                ds.x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}         
        ds.to_netcdf(path_run+f"{req[2][:-4]}.nc")


def mat2nc_mean_2D(path_run,path_sc="",run_file="run.sws"):
    """
    Save nc with dataset mean 2D (one nc for all output)
    Iterates over output requests and create single nc file with all mean quantities.
    Note: Only tested for multiple requests with one grid.
    path_sc is the path of SWASH source code if using developer mode.
    """
    import scipy.io as sio
    import xarray as xr
    _,reqn,_,_,_=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    out={}
    [out.update(sio.loadmat(path_run+mat[2])) for index,mat in enumerate(reqn) if reqn[index][1]=="mean" or 'XP' in reqn[index][3] or 'YP' in reqn[index][3]] # read only 'man' data
    out= dict(filter(lambda item:item[0] not in ('__header__', '__version__', '__globals__'),out.items())) # basic filter
    out= dict(filter(lambda item:not item[0][-1].isdigit(),out.items())) # filtering to remove 3D and instantaneous data (last character is digit)    
    ds=xr.Dataset(
                data_vars={key: (("y", "x"), value[::-1,:],{"standard_name": key,
                                                            "long_name": swash_dict.get(key,['']*2)[1],
                                                            "units": swash_dict.get(key,['']*3)[2]})
                        for key,value in out.items() if key not in ('Xp','Yp')},
                        attrs=attr_dict
                )

    if 'Yp' in out:
        ds=ds.assign_coords(y=out['Yp'][:,0][::-1]) #include cgrid here later on
        ds.y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
    if 'Xp' in out:
        ds=ds.assign_coords(x=out['Xp'][0,:])
        ds.x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}    
    print(f"Saving mean_2D.nc")
    ds.to_netcdf(path_run+f"mean_2D.nc")

def mat2nc_mean_3D(path_run,path_sc="",run_file="run.sws"):
    """
    Save nc with dataset mean 3D (one nc for all output)
    Iterates over output requests and create single nc file with all mean quantities.
    Note: Only tested for multiple requests with one grid (e.g. COMPGRID).
    path_sc is the path of SWASH source code if using developer mode.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    frame,reqn,nV,_,_=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    out={}
    [out.update(sio.loadmat(path_run+mat[2])) for index,mat in enumerate(reqn) if reqn[index][1]=="mean" or 'XP' in reqn[index][3] or 'YP' in reqn[index][3]] # read only 'mean' data
    out= dict(filter(lambda item:item[0] not in ('__header__', '__version__', '__globals__'),out.items())) # basic filter
    okv = list(set([i[:-1-int(np.log10(nV))] for i in dict(filter(lambda item:item[0][-1].isdigit(),out.items())).keys()])) # filter 3D data (ends with digit), will be datarray name
    okv_k=[[int(i[len(j):]) for i in out.keys() if i[:len(j)]==j] for j in okv] # will be k axis
    out_k= [[value for i,value in out.items() if i[:len(j)]==j] for j in okv] #will be datarray
    if out_k:
        ds=xr.Dataset(
                    data_vars={okv[index]: (("y", "x",["kc","ke"][min(okv_k[index])==0]), np.stack([v for v in value],axis=-1)[::-1],{"standard_name": okv[index],
                                                                "long_name": swash_dict.get(okv[index],['']*2)[1],
                                                                "units": swash_dict.get(okv[index],['']*3)[2]})
                            for index,value in enumerate(out_k)},
                            attrs=attr_dict
                        )
        if all([0 in okv_k[i] for i in range(len(okv_k))]): ds=ds.assign_coords(kc=okv_k[np.argmin([len(i) for i in okv_k])])
        if any([0 in okv_k[i] for i in range(len(okv_k))]): ds=ds.assign_coords(ke=list(filter(lambda i:0 in i,okv_k))[0])
        if 'Yp' in out:
            ds=ds.assign_coords(y=out['Yp'][:,0][::-1]) #include cgrid here later on
            ds.y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
        if 'Xp' in out:
            ds=ds.assign_coords(x=out['Xp'][0,:])
            ds.x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}   
        if 'kc' in ds:
            ds.kc.attrs = {"standard_name": 'kc',"long_name": 'Vertical axis (cell centre)', "units": '',"axis":"kc"}
        if 'ke' in ds:
            ds.ke.attrs = {"standard_name": 'ke',"long_name": 'Vertical axis (cell edge)', "units": '',"axis":"ke"}    
        print(f"Saving mean_3D.nc")
        ds.to_netcdf(path_run+f"mean_3D.nc")     

def mat2nc_ins_2D(path_run,path_sc="",run_file="run.sws"):
    """
    Save nc with dataset instantaneous 2D (one nc for each variable)
    Iterates over output requests and create one nc file for every variable.
    path_sc is the path of SWASH source code if using developer mode.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    frame,reqn,_,_,_=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    ivl=[] # var, y,  x and output (.mat)
    [[ivl.append([j,frame.get(i[0],['']*2)[0],frame.get(i[0],['']*2)[1],i[2]]) for j in  i[3] if j not in ['XP','YP','BOTL'] ] for i in reqn if i[1]=="inst"] # instantaneous variables list
    ivl = list(filter(lambda i:i[0][-1]!="K" and i[0]!="VZ",ivl)) #filter only 2D data, later on improve filter
    for i in ivl:
        out={}
        [out.update(sio.loadmat(path_run+mat[2])) for index,mat in enumerate(reqn) if i[3] == reqn[index][2]]
        out= dict(filter(lambda item:item[0] not in ('__header__', '__version__', '__globals__'),out.items())) # basic filter
        out= dict(filter(lambda item:swash_dict[item[0].split("_")[0]][0] in i[0],out.items())) # filter only relevant keys
        var = list(out.keys())[0].split("_")[0] #name of variable, should be unique
        ds=xr.Dataset(
                    data_vars={var: (("y", "x","t"), (np.stack(list(out.values()), axis=-1))[::-1],
                                                        {"standard_name": var,
                                                        "long_name": swash_dict[var][1],
                                                        "units": swash_dict[var][2]}),
                            },
                    coords={"y": i[1],
                            "x": i[2],
                            "t": [float(key[-10:-8])*60*60 +float(key[-8:-6])*60 +float(key[-6:-4]) + float(key[-3:])/1000 for key in out.keys()]
                        },
                        attrs=attr_dict
                    )
        ds.y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
        ds.x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}
        ds.t.attrs = {"standard_name": swash_dict['Tsec'][0],"long_name": swash_dict['Tsec'][1], "units": swash_dict['Tsec'][2],"axis":"t"}
        print(f"Saving ins_2D_{i[0]}_{(i[3].split('.'))[0].split('/')[-1]}.nc")
        ds.to_netcdf(path_run+f"ins_2D_{i[0]}_{i[3].split('.')[0].split('/')[-1]}.nc")


def mat2nc_ins_3D(path_run,path_sc="",run_file="run.sws"):
    """
    Save nc with dataset instantaneous 3D (one nc for each variable)
    Iterates over output requests and create one nc file for every variable.
    path_sc is the path of SWASH source code if using developer mode.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    frame,reqn,_,_,_=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    ivl=[] # var, y,  x and output (.mat)
    [[ivl.append([j,frame.get(i[0],['']*2)[0],frame.get(i[0],['']*2)[1],i[2]]) for j in  i[3] if j not in ['XP','YP','BOTL'] ] for i in reqn if i[1]=="inst"] # instantaneous variables list
    ivl = list(filter(lambda i:i[0][-1]=="K" or i[0]=="VZ",ivl)) #filter only 3D data, later on improve filter
    for i in ivl:
        out={}
        [out.update(sio.loadmat(path_run+mat[2])) for index,mat in enumerate(reqn) if i[3] == reqn[index][2]]
        out= dict(filter(lambda item:item[0] not in ('__header__', '__version__', '__globals__'),out.items())) # basic filter to exclude headings
        if i[0]!="VZ" and i[0]!="VELK":
            # only relevant keys for given request (works only for variables with Var_K01_time1_time2, thus not VZ or VELK or others to be found)            
            out= dict(filter(lambda item:len(item[0].split("_"))==4 and swash_dict["_".join("".join(filter(lambda i:not i.isdigit(),"_".join(item[0].split("_")[:2]))).split("_")[:2])][0]==i[0],out.items()))
            var = [j for i,j in enumerate(out.keys()) if i==0][0].split("_")[0]+"_k" #name of variable
            okv_k=list(set([int(i.split("_")[1][1:]) for i in out.keys()])) #will be k axis
            out_k=[dict(filter(lambda i:int(i[0].split("_")[1][1:])==k,out.items())) for k in okv_k] # filter by layer k
        elif i[0]=="VZ":
            out= dict(filter(lambda item:item[0][0]=="w",out.items())) # basic filter for vertical velocity
            var = "w" #name of variable
            okv_k=list(set([int(i.split("_")[0][1:]) for i in out.keys()])) #will be k axis
            out_k=[dict(filter(lambda i:int(i[0].split("_")[0][1:])==k,out.items())) for k in okv_k] # filter by layer k
        elif i[0]=="VELK": # not yet developed, code below is only for reference, e.g., requires separating x and y components before stacking
            out= dict(filter(lambda item:item[0][:3]=="Vel",out.items())) # basic filter for velocity vector
            var = "Vel_k" #name of variable
            okv_k=list(set([int(i.split("_")[2][1:]) for i in out.keys()])) #will be k axis
            out_k=[dict(filter(lambda i:int(i[0].split("_")[2][1:])==k,out.items())) for k in okv_k] # filter by layer k
        t = [float(key[-10:-8])*60*60 +float(key[-8:-6])*60 +float(key[-6:-4]) + float(key[-3:])/1000 for key in out_k[0].keys()]
        out_k=[np.stack(list(i.values()), axis=-1)[::-1] for i in out_k] # time stack
        out_k= np.stack([v for v in out_k],axis=-2) #layer stack
        ds=xr.Dataset(
                    data_vars={var: (("y", "x",["kc","ke"][min(okv_k)==0],"t"), out_k,
                                                        {"standard_name": var,
                                                        "long_name": swash_dict[var][1],
                                                        "units": swash_dict[var][2]}),
                            },
                    coords={"y": i[1],
                            "x": i[2],
                            ["kc","ke"][min(okv_k)==0] : okv_k,
                            "t": t
                            },
                    attrs=attr_dict
                    )
        ds.y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
        ds.x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}
        if 'kc' in ds:
            ds.kc.attrs = {"standard_name": 'kc',"long_name": 'Vertical axis (cell centre)', "units": '',"axis":"kc"}
        if 'ke' in ds:
            ds.ke.attrs = {"standard_name": 'ke',"long_name": 'Vertical axis (cell edge)', "units": '',"axis":"ke"}
        ds.t.attrs = {"standard_name": swash_dict['Tsec'][0],"long_name": swash_dict['Tsec'][1], "units": swash_dict['Tsec'][2],"axis":"t"}        
        print(f"Saving ins_3D_{i[0]}_{i[3].split('.')[0].split('/')[-1]}.nc")
        ds.to_netcdf(path_run+f"ins_3D_{i[0]}_{i[3].split('.')[0].split('/')[-1]}.nc")

def mat2nc_ins_table(path_run,path_sc="",run_file="run.sws"):
    """
    Save nc with dataset instantaneous table (does not work for vector output)
    It saves one netcdf file per line request.
    """
    import scipy.io as sio
    import xarray as xr
    import numpy as np
    _,_,_,reqtable,point=load_req(path_run,run_file=run_file)
    swash_dict,attr_dict=swashdict(path_sc=path_sc)
    swash_dict_rev={j[0]:[i]+j[1:] for i,j in swash_dict.items()}
    for tab in reqtable: #name, of grid type (mean or inst), name of tbl, list of requested variables
        out={}
        if  open(path_run+tab[2]).readlines()[0][:5]!="SWASH": # no header 
            for idx,var in enumerate(tab[3]):
                out[var]= np.transpose(np.reshape(np.array([float(j[idx]) for j in [i.strip().split()
                    for i in open(path_run+tab[2]).readlines()]]),
                    (len(open(path_run+tab[2]).readlines())//
                    len(point[tab[0]][0]),len(point[tab[0]][0])))) # station,time
        else: # with header --> reads number of dimensions in vertical for each variable
            ndim=[int(i) for i in [i.split()[0] for i in open(path_run+tab[2]).readlines()][8:8+4*(len(tab[-1])-1)+1][::4]]
            outarray=np.concatenate([[float(k) for k in j.strip().split()] for j in open(path_run+tab[2]).readlines()[8+4*(len(tab[-1])-1)+2:]]) # one line
            outlist=[]
            for idx,var in enumerate(tab[3]): # parse into text file and append list for each variable
                outlist.append([])
                i=0
                for j in range(idx): i+=ndim[j]
                while i<len(outarray):
                    for j in range(ndim[idx]):
                        outlist[idx].append(outarray[i])
                        i+=1
                    i+=int(np.dot(ndim,np.ones(len(ndim)))-ndim[idx])
                if ndim[idx]==1:
                    out[var]=np.transpose(np.reshape(np.array(outlist[idx]),(len(outlist[idx])//len(point[tab[0]][0]),len(point[tab[0]][0]))))
                else:
                    out[var]=np.transpose(np.reshape(np.array(outlist[idx]),(len(outlist[idx])//len(point[tab[0]][0])//ndim[idx],len(point[tab[0]][0]),ndim[idx])),(1,0,2))
        
        var_2D = [i for i,j in zip(tab[-1],ndim) if j==1]
        ds_2D=xr.Dataset(
                    data_vars={swash_dict_rev.get(key,[''])[0]: (("stations", "time"), value,{"standard_name": swash_dict_rev.get(key,[''])[0],
                                                                "long_name": swash_dict_rev.get(key,['']*2)[1],
                                                                "units": swash_dict_rev.get(key,['']*3)[2]})
                            for key,value in out.items() if key in var_2D},
                            attrs=attr_dict
                    )
        ds_3D=xr.Dataset(
                    data_vars={swash_dict_rev.get(key,[''])[0]: (("stations","time","k"), value,{"standard_name": swash_dict_rev.get(key,[''])[0],
                                                                "long_name": swash_dict_rev.get(key,['']*2)[1],
                                                                "units": swash_dict_rev.get(key,['']*3)[2]})
                            for key,value in out.items() if key not in ('XP','YP','TIME','TSEC')
                                                            and key not in var_2D},
                            attrs=attr_dict
                    )
        ds=xr.merge([ds_2D,ds_3D])
        if 'TIME' in out.keys():
            t="Time"
            ds=ds.assign_coords(time=out["TIME"][0,:])
            ds.time.attrs = {"standard_name": swash_dict[t][0],"long_name": swash_dict[t][1], "units": swash_dict[t][2],"axis":"time"}
        elif 'TSEC' in out:
            t="Tsec"
            ds=ds.assign_coords(time=out["TSEC"][0,:])
            ds.time.attrs = {"standard_name": swash_dict[t][0],"long_name": swash_dict[t][1], "units": swash_dict[t][2],"axis":"time"}
        ds['station_y']=(("stations"),point[tab[0]][0])
        ds['station_x']=(("stations"),point[tab[0]][1])
        ds.station_y.attrs = {"standard_name": swash_dict['Yp'][0],"long_name": swash_dict['Yp'][1], "units": swash_dict['Yp'][2],"axis":"Y"}
        ds.station_x.attrs = {"standard_name": swash_dict['Xp'][0],"long_name": swash_dict['Xp'][1], "units": swash_dict['Xp'][2],"axis":"X"}    
        print(f"Saving ins 2D table {tab[2].split('.')[0]}.nc")
        ds.to_netcdf(path_run+f"{tab[2].split('.')[0]}.nc")

if __name__ == '__main__':
    pass