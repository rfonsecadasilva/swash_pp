"""
Created on May 04 2021
Write nc files from SWASH gridded mat output runs.
@author: rfonsecadasilva
"""

def swashdict(path_sc=""):
    """
    Return dictionary with SWASH output metadata from source code.

    Args:
        path_sc (str, optional): path of SWASH source code if using developer mode (otherwise it parses SWASH 8.01 online version). Defaults to "".
    
    Returns:
        swash_dict (dict): dictionary with SWASH output metadata (keyword, name, and units; read from SwashInit.ftn90)
        attr_dict (dict): dictionary with SWASH version metadata
    """
    def swashunit(init):
        """
        Return dict swash_units (e.g., 'uh' corresponds to m)
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
        tmpfile.write(urllib.request.urlopen("https://swash.sourceforge.io/download/zip/swash-8.01.tar.gz").read())
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
    Parse SWASH input file and return relevant information for creating nc files.

    Args:
        path_run (str): path where run file (*.sws) is located.
        run_file (str, optional): run file name. Defaults to "run.sws".
    
    Returns:
        frame (list): list with grids (requested with BLOCK command) - id, y, and x
        reqn (list): list with gridded requests - id, type (mean or inst), and list of requested variables
        nv (int): number of vertical layers
        reqtable (list): list with table  requests - id, type (mean or inst), and list of requested variables
        frame (list): list with table points (requested with POI command) - id, y, and x
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
    if cgrid: frame={**frame,**cgrid}
    return frame,reqn,nV,reqtable,point

def mat2nc_all(path_run,path_sc="",run_file="run.sws"):
    """
    Convert all gridded matlab output into nc files

    Args:
        path_run (str): path where run file (*.sws) is located.
        path_sc (str, optional): path of SWASH source code if using developer mode (otherwise it parses SWASH 8.01 online version). Defaults to "".
        run_file (str, optional): run file name. Defaults to "run.sws".  
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