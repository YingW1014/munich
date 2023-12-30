#! /usr/bin/env python

import math, re, datetime, sys
import numpy as np
from io_tools import *
import netCDF4

# Read background binary files from Polair3D output
def read_bkgd_bin(bkgd_species, indir, Nt, Nx, Ny, Nz):
    # Put data in a dict
    bkgd_data = {}
    for species in bkgd_species:
        infile = indir + species + '.bin'
        bkgd_data[species] = np.memmap(infile, dtype='float32', mode='r',
                                       shape=(Nt, Nz, Ny, Nx))[:, 0, :, :]

    return bkgd_data


# Set background for streets
def set_bkgd_bin(street_list, bkgd_data, current_date, date_min, delta_t,
                  Nt, x_min, y_min, delta_x, delta_y, Nx, Ny):
    # Get index of current date
    c_id = int((current_date - date_min).total_seconds() / delta_t)
    if c_id < 0 or c_id >= Nt:
        sys.exit('ERROR: background data not available for this date.')

    for i in range(len(street_list)):
        # Get cell indices (X, Y) in the Polair3D grid
        # Xid = (street_list[i].lon_cen - x_min) / delta_x
        # Yid = (street_list[i].lat_cen - y_min) / delta_y
        # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
        # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)
        Xid, Yid = get_polair_id(street_list[i].lon_cen, street_list[i].lat_cen,
                                 x_min, y_min, delta_x, delta_y, Nx, Ny)
        for key, value in bkgd_data.items():
            street_list[i].bkgd[key] = value[c_id, Yid, Xid]

# Append background to binary files
def append_bkgd_data(street_list, bkgd_files):
    for i, var in enumerate(street_list[0].bkgd.keys()):
        data = np.zeros((len(street_list)), 'float')
        for j in range(len(street_list)):
            data[j] = street_list[j].bkgd[var]
        append_binary(data, bkgd_files[var])


# Read the background concentration from a text file.
def get_background_concentration(input_file, current_date, street_list):


    hasBackground = False
    input_background = open(input_file)
    header = input_background.readline()
    nstreet = 0
    for line in input_background.readlines():
        line = line.replace('\n','')
        line_info = [x for x in re.split('\t| ', line) if x.strip() != '']
        str_times = line_info[0]
        year = int(str_times[0:4])
        month = int(str_times[5:7])
        day = int(str_times[8:10])
        hour = int(str_times[11:13])

        background_date = datetime.datetime(year,month,day,hour)
        if current_date == background_date:
            print("Background concentration data are available for the date ", current_date)
            o3 = float(line_info[1])
            no2 = float(line_info[2])
            no = float(line_info[3])
            hasBackground = True
            break

    if (hasBackground == False):
        print("Error: background concentration data are not available.")
    else:   
        for s in range(len(street_list)):
            street = street_list[s]
            street.bkgd['O3'] = o3
            street.bkgd['NO2'] = no2
            street.bkgd['NO'] = no

    return 0
        
def get_polair_id(lon, lat, x_min, y_min, dx, dy, Nx, Ny):
    # Get cell indices (X, Y) in the Polair3D grid
    # Xid = (street_list[i].lon_cen - x_min) / delta_x
    # Yid = (street_list[i].lat_cen - y_min) / delta_y
    # Xid = int(Xid)-1 if Xid%1 < 0.5 else int(Xid)
    # Yid = int(Yid)-1 if Yid%1 < 0.5 else int(Yid)

    Xid = max(int((lon - x_min + dx / 2.) / dx), 0)
    Xid = min(Xid, Nx - 1)
    Yid = max(int((lat - y_min + dy / 2.) / dy), 0)
    Yid = min(Yid, Ny - 1)

    return Xid, Yid

def get_polair_ind(polair_lon, polair_lat, street):
    min_distance = 99.
    i = 0
    j = 0
    for lon in polair_lon:
        for lat in polair_lat:
            distance = np.sqrt(pow((lon - street.lon_cen), 2.0) + 
                               pow((lat - street.lat_cen), 2.0))
            if distance < min_distance:
                indx = i
                indy = j
                min_distance = distance
            j = j + 1
        j = 0    
        i = i + 1
    return indx, indy

def get_chimere_background_concentration(current_date, street_list, melchior_spec_list,molar_mass_melchior2,chimout_dir,chimout_lab) :

    import netCDF4
    import os,sys

    str_date = current_date.strftime("%Y%m%d")
    input_file=chimout_dir+'/out.'+str_date+'00_'+chimout_lab+'.nc'
    if not os.path.isfile(input_file) :
       print(('CHIMERE background conditions are requested but the file is not found: '+str(input_file)))
       sys.exit()

    nc = netCDF4.Dataset(input_file, 'r')
    chim_times = nc.variables["Times"][:]
    lons = nc.variables["lon"][:]
    lats = nc.variables["lat"][:]
    new_spec_list=[] 
    # New melchior species list from what is actually present in the CHIMERE out file
    for spec in melchior_spec_list:
        if spec in list(nc.variables.keys()):
           new_spec_list.append(spec)
        else:
           print(('Warning!! '+spec+' not found in CHIMERE output file'))

    # Transform CHIMERE date-time to python datetime

    N=chim_times.shape[0] #number of hours to parse
    times=[]
    for i in range(N): #loop over hours
            s1,s2,s3,s4=chim_times[i,0:4]
            YEAR=s1+s2+s3+s4
            s1,s2=chim_times[i,5:7]
            MONTH=s1+s2
            s1,s2=chim_times[i,8:10]
            DAY=s1+s2
            s1,s2=chim_times[i,11:13]
            HOUR=s1+s2
            times.append(datetime.datetime(int(YEAR),int(MONTH),int(DAY),int(HOUR)))
    times=np.array(times)

    for spec in new_spec_list:
        for s in range(len(street_list)):
            street = street_list[s]
            street.background[spec]=0.0


    hasBackground = False
    for t in range(len(times)):
        background_date = times[t]
        if current_date == background_date:
            print("Background data are available for the date ", current_date)
            ind_t = t
            hasBackground = True

    if (hasBackground == False):
        print("Error: background data are not available")
        sys.exit()

    print((lons.shape))
    nx = lons.shape[1]
    ny = lons.shape[0]

    for s in range(len(street_list)):
       street = street_list[s]
       lat1 = street.lat_cen
       lon1 = street.lon_cen
       init_length = 9999.0
       for i in range(nx):
           for j in range(ny):
               lat2 = lats[j, i]
               lon2 = lons[j, i]
               length = distance_on_unit_sphere(lat1, lon1, lat2, lon2)
               if length < init_length:
                   init_length = length
                   ind_i = i
                   ind_j = j
       print((ind_i,ind_j))
       tem2=nc.variables['tem2'][ind_t,ind_j,ind_i] #Kelvin
       psfc=nc.variables['pres'][ind_t,0,ind_j,ind_i] #Pascal

       if psfc > 0 :
          molecular_volume=22.41 * (tem2/273.)  * (1013*10**2)/psfc
       else :
          molecular_volume=22.41

       for spec in new_spec_list:
           conc_ppb=nc.variables[spec][ind_t,0,ind_j,ind_i]
           molecular_mass=molar_mass_melchior2[spec]
           ppb2ug=molecular_mass / molecular_volume   #
           if spec=='NO2' : print((s,ind_t,0,ind_j,ind_i,conc_ppb,conc_ppb * ppb2ug))
           street.background[spec]=conc_ppb * ppb2ug

    return 0     #street.background in ug/m3


def get_wrfchem_background_concentration(current_date, street_list, wrfchem_gas_spec_list, \
                    molar_mass_melchior2, input_dir,  wrf_config) : # Developed by YingWang (YW)

    import math
    
    content =  [("file_type", "[type]", "Int"), \
                ("wrfout_prefix", "[type]", "String"), \
                ("filename","[type]","String"), \
                ('time_step_wrf', '[type]', 'Float'), \
                ('option_lmo', '[option]', 'String'), \
                ('option_surface_wind', '[option]', 'String')
    ]

    print("WRF configuration file: ", wrf_config)
    config = talos.Config(wrf_config, content)
    wrfout_prefix = config.wrfout_prefix

#    from scipy.io import netcdf

    str_date = current_date.strftime("%Y-%m-%d")
    input_file = input_dir + wrfout_prefix + "_" + str_date + "-00:00:00"
    if not os.path.isfile(input_file) :
       print(('WRF/Chem background conditions are requested but the file is not found: '+str(input_file)))
       sys.exit()

#    nc = netcdf.netcdf_file(input_file, 'r')
    nc = netCDF4.Dataset(input_file, 'r')  # By YW
    wrf_times = nc.variables["Times"][:]
    lons = nc.variables["XLONG"][:]
    lats = nc.variables["XLAT"][:]
    ALT = nc.variables["ALT"][:] # inverse density
    
    wrf_aero_var_p1 = ['so4','nh4','no3','na','cl']
    wrf_aero_var_p2 = ['ec','p25','orgpa']
    MUNICH_aero_var_p1 = ['PSO4_','PNH4_','PNO3_','PNA_','PHCL_']
    MUNICH_aero_var_p2 = ['PBC_','P25_','PPOAlP_']
    MUNICH_soa_var = ['PSOAlP_','PSOAmP_','PSOAhP_']
    
    # Redistribute mass concentrations from lognormal distribution to size bin
    Bin_bounds = [0.01,0.0398,0.1585,0.4,1.0,2.5115,10.]
    
    # Standard deviation and volume (mass) geometric mean diameter in um, 
    # Parameters are collected from chem/module_data_sorgam_vbs.F ;from (Ackermann et al., Atm. Env., 1998)
    
    sg_i = 1.70 
    dgv_i = 0.03 
    dgn_i = 0.01
    sg_j = 2.0 
    dgv_j = 0.3
    dgn_j = 0.07
    sg_c = 2.2
    dgv_c = 6.0
    dgn_c = 1.0

    new_spec_list=[] 

    # New WRF/Chem gas species list from what is actually present in the WRF/Chem out file
    for spec in wrfchem_gas_spec_list:
        if spec.lower() in list(nc.variables.keys()):
           new_spec_list.append(spec)
        else:
           print(('Warning!! '+spec+' not found in WRF/Chem output file'))

    # Transform WRF/Chem date-time to python datetime
    N=wrf_times.shape[0] #number of hours to parse
    times=[]
    for i in range(N): #loop over hours
            s1,s2,s3,s4=wrf_times[i,0:4]
            YEAR=s1+s2+s3+s4
            s1,s2=wrf_times[i,5:7]
            MONTH=s1+s2
            s1,s2=wrf_times[i,8:10]
            DAY=s1+s2
            s1,s2=wrf_times[i,11:13]
            HOUR=s1+s2
            times.append(datetime.datetime(int(YEAR),int(MONTH),int(DAY),int(HOUR)))
    times=np.array(times)

    for spec in new_spec_list:
        for s in range(len(street_list)):
            street = street_list[s]
            street.background[spec]=0.0


    hasBackground = False
    for t in range(len(times)):
        background_date = times[t]
        if current_date == background_date:
            print("WRF/Chem background data are available for the date ", current_date)
            ind_t = t
            hasBackground = True

    if (hasBackground == False):
        print("Error: WRF/Chem background data are not available")
        sys.exit()

    #print((lons.shape))
    nx = lons.shape[2]
    ny = lons.shape[1]
    
    # Calculate background concentrations from WRF/Chem
    for s in range(len(street_list)):
        street = street_list[s]
        lat1 = street.lat_cen
        lon1 = street.lon_cen
        init_length = 9999.0
        for i in range(nx):
            for j in range(ny):
                lat2 = lats[0, j, i]
                lon2 = lons[0, j, i]
                length = distance_on_unit_sphere(lat1, lon1, lat2, lon2)
                if length < init_length:
                    init_length = length
                    ind_i = i
                    ind_j = j
        #print((ind_i,ind_j))
        tem2=nc.variables['T2'][ind_t,ind_j,ind_i] #Kelvin
        psfc=nc.variables['PSFC'][ind_t,ind_j,ind_i] #Pascal

        if psfc > 0 :
            molecular_volume=22.41 * (tem2/273.)  *  1.013*math.pow(10,5)/psfc
        else :
            molecular_volume=22.41

        # Calculate background concentrations for gaseous species
        for spec in new_spec_list:
            conc_ppm=nc.variables[spec.lower()][ind_t,0,ind_j,ind_i]
            molecular_mass=molar_mass_melchior2[spec]
            ppm2ug= 1000 * molecular_mass / molecular_volume   #
            # if spec=='NO2' : print((s,ind_t,0,ind_j,ind_i,conc_ppb,conc_ppb * ppb2ug))
            street.background[spec]=conc_ppm * ppm2ug
            # print('Unit conversion for ' + spec + ' from ' + str(conc_ppm) + ' to ' + str(street.background[spec]))
        print('Complete background for gases')

        # Calculate background concentrations for aerosol species
        # Group1 
        for k in range(len(MUNICH_aero_var_p1)):
            MUNICH_species = MUNICH_aero_var_p1[k]
            WRF_species = wrf_aero_var_p1[k]
            for i in range(6):
                name_MUNICH = MUNICH_species + str(i)
                #print(name_MUNICH)
                delta_d = Bin_bounds[i+1]-Bin_bounds[i]
                d = (Bin_bounds[i+1] + Bin_bounds[i])/2
                
                # In aitke mode
                MLi = (nc.variables[WRF_species + 'ai'][ind_t,0,ind_j,ind_i] + nc.variables[WRF_species + 'cwi'][ind_t,0,ind_j,ind_i]) / ALT[ind_t,0,ind_j,ind_i]
                Eq_part1_i = MLi * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_i))
                Eq_part2_i = math.exp(-math.pow(math.log(d/dgv_i),2) / (2*math.pow(math.log(sg_i),2)))
                # In accumation mode
                MLj = (nc.variables[WRF_species + 'aj'][ind_t,0,ind_j,ind_i] + nc.variables[WRF_species + 'cwj'][ind_t,0,ind_j,ind_i]) / ALT[ind_t,0,ind_j,ind_i]
                Eq_part1_j = MLj * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_j))
                Eq_part2_j = math.exp(-math.pow(math.log(d/dgv_j),2) / (2*math.pow(math.log(sg_j),2)))
                
                street.background[name_MUNICH] = Eq_part1_i * Eq_part2_i + Eq_part1_j * Eq_part2_j


        for k in range(len(MUNICH_aero_var_p2)):
            MUNICH_species = MUNICH_aero_var_p2[k]
            WRF_species = wrf_aero_var_p2[k]
            for i in range(6):
                name_MUNICH = MUNICH_species + str(i)
                #print(name_MUNICH)
                delta_d = Bin_bounds[i+1]-Bin_bounds[i]
                d = (Bin_bounds[i+1] + Bin_bounds[i])/2

                # In aitke mode
                MLi = (nc.variables[WRF_species + 'i'][ind_t,0,ind_j,ind_i] + nc.variables[WRF_species + 'cwi'][ind_t,0,ind_j,ind_i]) / ALT[ind_t,0,ind_j,ind_i]
                Eq_part1_i = MLi * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_i))
                Eq_part2_i = math.exp(-math.pow(math.log(d/dgv_i),2) / (2*math.pow(math.log(sg_i),2)))
                # In accumation mode
                MLj = (nc.variables[WRF_species + 'j'][ind_t,0,ind_j,ind_i] + nc.variables[WRF_species + 'cwj'][ind_t,0,ind_j,ind_i]) / ALT[ind_t,0,ind_j,ind_i]
                Eq_part1_j = MLj * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_j))
                Eq_part2_j = math.exp(-math.pow(math.log(d/dgv_j),2) / (2*math.pow(math.log(sg_j),2)))

                street.background[name_MUNICH] = Eq_part1_i * Eq_part2_i + Eq_part1_j * Eq_part2_j


        for k in range(len(MUNICH_soa_var)):
            MUNICH_species = MUNICH_soa_var[k]
            if k == 0:
                temp1 = nc.variables['asoa1i'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa1i'][ind_t,0,ind_j,ind_i]
                temp2 = nc.variables['asoa1j'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa1j'][ind_t,0,ind_j,ind_i]
            elif k==1:
                temp1 = nc.variables['asoa2i'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa2i'][ind_t,0,ind_j,ind_i] + nc.variables['asoa3i'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa3i'][ind_t,0,ind_j,ind_i]
                temp2 = nc.variables['asoa2j'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa2j'][ind_t,0,ind_j,ind_i] + nc.variables['asoa3j'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa3j'][ind_t,0,ind_j,ind_i]
            elif k==2:
                temp1 = nc.variables['asoa4i'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa4i'][ind_t,0,ind_j,ind_i]
                temp2 = nc.variables['asoa4j'][ind_t,0,ind_j,ind_i] + nc.variables['bsoa4j'][ind_t,0,ind_j,ind_i]
            
            for i in range(6):
                name_MUNICH = MUNICH_species + str(i)
                #print(name_MUNICH)
                delta_d = Bin_bounds[i+1]-Bin_bounds[i]
                d = (Bin_bounds[i+1] + Bin_bounds[i])/2
                
                # In aitke mode
                Eq_part1_i = temp1 / ALT[ind_t,0,ind_j,ind_i] * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_i))
                Eq_part2_i = math.exp(-math.pow(math.log(d/dgv_i),2) / (2*math.pow(math.log(sg_i),2)))

                # In accumation mode
                Eq_part1_j = temp2 / ALT[ind_t,0,ind_j,ind_i] * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_j))
                Eq_part2_j = math.exp(-math.pow(math.log(d/dgv_j),2) / (2*math.pow(math.log(sg_j),2)))

                street.background[name_MUNICH] = Eq_part1_i * Eq_part2_i + Eq_part1_j * Eq_part2_j
                
        for i in range(6):
            name_MUNICH = 'PMD_' + str(i)
            
            delta_d = Bin_bounds[i+1]-Bin_bounds[i]
            d = (Bin_bounds[i+1] + Bin_bounds[i])/2
            
            # In coarse mode
            MLc = (nc.variables['soila'][ind_t,0,ind_j,ind_i] + nc.variables['antha'][ind_t,0,ind_j,ind_i] + nc.variables['seas'][ind_t,0,ind_j,ind_i] + 
                   nc.variables['soilcw'][ind_t,0,ind_j,ind_i] + nc.variables['anthcw'][ind_t,0,ind_j,ind_i] + nc.variables['seascw'][ind_t,0,ind_j,ind_i]) / ALT[0,0,0,0]
            Eq_part1_c =  MLc * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_c))
            Eq_part2_c = math.exp(-math.pow(math.log(d/dgv_c),2) / (2*math.pow(math.log(sg_c),2)))
            
            street.background[name_MUNICH] = Eq_part1_c * Eq_part2_c
        
        print('Complete background for aerosols')
                    
        # Calculate number concentrations
        for i in range(6):
            name_MUNICH = 'Number_' + str(i)

            delta_d = Bin_bounds[i+1]-Bin_bounds[i]
            d = (Bin_bounds[i+1] + Bin_bounds[i])/2
            
            NLi = nc.variables['nu0'][ind_t,0,ind_j,ind_i] + nc.variables['nu0cw'][ind_t,0,ind_j,ind_i]
            NLj = nc.variables['ac0'][ind_t,0,ind_j,ind_i] + nc.variables['ac0cw'][ind_t,0,ind_j,ind_i]
            NLc = nc.variables['corn'][ind_t,0,ind_j,ind_i] + nc.variables['corncw'][ind_t,0,ind_j,ind_i]
            
            # In aitke mode
            Eq_part1_i =  NLi * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_i))
            Eq_part2_i = math.exp(-math.pow(math.log(d/dgn_i),2) / (2*math.pow(math.log(sg_i),2)))
            # In accumation mode
            Eq_part1_j = NLj * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_j))
            Eq_part2_j = math.exp(-math.pow(math.log(d/dgn_j),2) / (2*math.pow(math.log(sg_j),2)))
                    
            Eq_part1_c =  NLc * delta_d / (d*math.sqrt(2*math.pi) * math.log(sg_c))
            Eq_part2_c = math.exp(-math.pow(math.log(d/dgn_c),2) / (2*math.pow(math.log(sg_c),2)))
        
            street.background[name_MUNICH] = Eq_part1_i * Eq_part2_i + Eq_part1_j * Eq_part2_j + Eq_part1_c * Eq_part2_c
        
        print('Complete background for number')

    return 0     #street.background in ug/m3