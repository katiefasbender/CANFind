#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# transform_to_mpc80.py is a script that reads fits files in a subdir (hgroup_<subdir#>) with CANFind output,
# grabs the tracklets, and records them in the Minor Planet Center's 80-column format,
# explained in no detail below (find more at https://www.minorplanetcenter.net/iau/info/OpticalObs.html),
# then writes the new tracklet table to a text file.  This is the file type that digest2 and MPC require.
# It also queries any necessary information? We'll see about that. The tracklet measurements are also
# written to a FITS file, to later be concatenated into the tracklet catalog by cf_tracklet_cat_creator.py

# command format:  $python /path/to/transform_to_mpc80.py <subdir#> 

# From the MPC (that website listed above):
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For Minor Planets:   (assume this for all tracklets, initially )
#    1 -  5       A5     Packed minor planet number
#    6 - 12       A7     Packed provisional designation, or a temporary designation
#    13           A1     Discovery asterisk
#   For Comets:
#    1 -  4       I4     Periodic comet number
#    5            A1     Letter indicating type of orbit
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#   For Natural Satellites:
#    1            A1     Planet identifier [Only if numbered]
#    2 -  4       I3     Satellite number  [Only if numbered]
#    5            A1     "S"
#    6 - 12       A7     Provisional or temporary designation
#    13           X      Not used, must be blank
#---------------------------
#   Columns     Format   Use
#---------------------------
#   For all: (columns 14-80)
#    14            A1    Note 1
#    15            A1    Note 2
#    16 - 32             Date of observation
#    33 - 44             Observed RA (J2000.0)
#    45 - 56             Observed Dec (J2000.0)
#    57 - 65       9X    Must be blank
#    66 - 71    F5.2,A1  Observed magnitude and band
#    72 - 77       X     Must be blank
#    78 - 80       A3    Observatory code

# Example Row:
#
#
#

#---------------------------------------------------------------------------------
# Imports
#---------------------------------------------------------------------------------
import matplotlib
import sys
import os
import numpy as np
from astropy.table import Table,Column,join,vstack
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.time import Time
from dl import queryClient as qc
import healpy as hp

#---------------------------------------------------------------------------------
# Functions
#---------------------------------------------------------------------------------
def query_fix(quer,dtps):
    bc=[]
    bdd=quer.splitlines()
    for ag in bdd:
        bc.append(ag.split(","))
    return(Table(np.array(bc[1:]),names=np.array(bc[0]),dtype=dtps))

def makedir(dir):
    '''Makes a directory with name "dir" if it does not exist.
    Arguments:
    ----------
    dir (str)
            directory name
    Returns:
    --------
    None; directory "dir" is created if not already there
    '''
    if not os.path.exists(dir):
        os.mkdir(dir)

def calc_mpc80(obs,scope_names,scope_codes):
    '''Calculates necessary info for MPC 80-column format line for observation.
    Arguments:
    ----------
    obs (astropy row)
            a tracklet observation
    scope_names (list)
            names of telescope/observatories used
    scope_codes (list)
            codes for "                          "
    Returns:
    --------
    None; directory "dir" is created if not already there
    '''
#Date of Observation: "YYYY MM DD.dddddd"
    tim=Time(float(obs['mjd']),format="mjd")
    yyyy=str(tim.datetime.year).zfill(4)
    mm=str(tim.datetime.month).zfill(2)
    dd='{:.6f}'.format(round((tim.datetime.day+(tim.datetime.hour+
        ((tim.datetime.minute+((tim.datetime.second+(tim.datetime.microsecond/1000000))/60))/60))/24),6))
    obs_date=" ".join([yyyy,mm,dd.zfill(9)])
#Observed RA: "HH MM SS.ddd",
    co=SkyCoord(obs['ra'],obs['dec'],frame="icrs",unit="degree")
    hhr=int(co.ra.hour)
    mmr=int((co.ra.hour-hhr)*60)
    ssr='{:.3f}'.format(round(((co.ra.hour-hhr)*60-mmr)*60,3))
    obs_ra=" ".join([str(hhr).zfill(2),str(mmr).zfill(2),str(ssr).zfill(6)])
#Observed Dec: "sDD MM SS.dd",
    dd=int(co.dec.deg)
    sine=str("+")
    if np.sign(dd)==-1:
        sine=str("-")
    mm=int((co.dec.deg-dd)*60)
    ss='{:.2f}'.format(abs(round((((co.dec.deg-dd)*60-mm)*60),2)))
    obs_dec="".join([sine,str(abs(dd))," ",str(abs(mm)).zfill(2)," ",str(ss).zfill(5)])
#Packed Provisional Number:
    #p1="A"
    #for cent in range(0,len(cent_numbers)):
    #   print(int(yyyy[:3]))
    #   if int(yyyy[:3])==cent_numbers[cent]:
    #       p1=cent_letters[cent]
    #p_des=str(pix).zfill(6) #my temporary designation number(the healpix number)
#Observed Magnitude and Band: mm.dd,I
    obs_mag="".join([str(round(obs['mag_auto'],2)).zfill(5),obs['filter'][0]])
#Observatory code:
    instrument=obs['instrument']
    for ins in range(0,len(scope_names)):
            if instrument==scope_names[ins]: 
                obs_code=scope_codes[ins]
    return(obs_date,obs_ra,obs_dec,obs_mag,obs_code)

def make_80col(obs_num,t_num,obs_date,obs_ra,obs_dec,obs_mag,obs_code):
    '''Formats observation info into MPC 80-column line for Find_Orb use.
    Arguments:
    ----------
    obs_num (int)
            a unique observation number
    t_num (int)
            a unique tracklet number
    the rest is output from calc_mpc80()
    Returns:
    --------
    None; directory "dir" is created if not already there
    '''
    c01t05=str(obs_num).zfill(6) #columns 1-5 observation number, left-padded with 0s
    c06t12=str(t_num).zfill(7)   #columns 6-12 packed provisional number (pix for now)
    c13t15=str("   ")            #columns 13-15 3X
    c16t32=str(obs_date)         #columns 16-32 Date of observation
    c33t44=str(obs_ra)           #columns 33-44 Observed RA (J2000.0)
    c45t56=str(obs_dec)          #columns 45-56 Observed Decl. (J2000.0)
    c57t65=str("         ")      #columns 57-65 9X     Must be blank
    c66t71=str(obs_mag)          #columns 66-71 F5.2,A1   Observed magnitude and band
    c72t77=str("      ")         #columns 72-77 6X      Must be blank
    c78t80=str(obs_code)         #columns 78-80 A3     Observatory code
    the_line="".join([c01t05,c06t12,c13t15,c16t32,c33t44,c45t56,c57t65,c66t71,c72t77,c78t80,str("\n")]) #the 80col line to write!
    return the_line

#---------------------------------------------------------------------------------
# Main Code
#---------------------------------------------------------------------------------
if __name__ == "__main__":

    # --Initiating & Defining--
    qc.set_timeout_request(1800)                # set datalab query timeout 
    basedir = "/home/x25h971/canfind_dr2/"      # base location of operations
    scope_codes=["V00","695","W84"]             # telescope information
    scope_names=["ksb","k4m","c4d"]
    cent_letters=["I","J","K"]                  # century of observation
    cent_numbers=[int(18),int(19),int(20)]      
    qnames=["measid","objectid","exposure"]     # nsc.meas columns to query if missing 
    qdtypes=['U','U','U']
    file_exposures=Table.read(basedir+"files/"+"exp_inst.fits",format="fits")  # nsc_dr1.exposure table (exp_inst.fits)

    # Set up table for hgroup nsc.meas info; this table will contain all tracklet obs info for HP in the hgroup subdir
    # also - these are the columns each HP should have.  
    ns=["mjd","ra","dec","measid","objectid","mag_auto","magauto_err","filter","raerr","decerr","exposure"]
    dts=["float64","float64","float64","U","U","float64","float64","U","float64","float64","U"]
    cat_mmts=Table(names=ns,dtype=dts)

    # Set up the text file for subdir to store MPC-formatted data (MPC80 file)
    subdir=sys.argv[1] #the subdirectory number containing the files to combine (hgroup_subdir)
    makedir(basedir+"concat_files/hgroup_%s"%subdir)
    file80=open("cf_tracklet_obs_%s.txt" % subdir,"a+") 
    
    # --Loop through HP files in subdir--
    # to query missing data, combine tables, & format tracklet mmt data into 80-col format
    t_num=0 #counter for tracklet_number (every tracklet in this subdir hgroup will have a unique number)
    for root,dirs,files in os.walk("hgroup_%s"%str(subdir)):
        for name in files:
            if os.stat(os.path.join(root,name)).st_size!=0 and len(fits.open(os.path.join(root,name)))>1:
                
                # read file, get HP number & tracklet info 
                hdul=fits.open(os.path.join(root,name))
                dat=Table(hdul[1].data)
                hdul.close()
                pix=(str(name).split("_")[-1]).split(".")[0]
                bool_out=dat['cluster_label']!=-1
                t_out=dat[bool_out] #just the tracklet mmts 

                # if there are tracklets in this hpfile,
                if len(t_out)>0 and (not (str(pix).zfill(6)) in file80.read()):

                # --Query missing columns from nsc.meas--
                    if "exposure" not in t_out.colnames:
                        RA=hp.pix2ang(128,int(pix),lonlat=True)[0]
                        DEC=hp.pix2ang(128,int(pix),lonlat=True)[1]
                        nbs=hp.get_all_neighbours(512,RA,DEC,lonlat=True) #8 nearest neighbors to the cooordinates for nside=512
                        coords=hp.pix2ang(512,nbs,lonlat=True) #center coordinates for the 8 nearest neighbors
                        fpix=np.unique(hp.ang2pix(256,coords[0],coords[1],lonlat=True))
                        qtext="".join(["SELECT meas.measid,meas.objectid,meas.exposure FROM nsc_dr1.meas as meas JOIN nsc_dr1.object as obj on obj.id=meas.objectid WHERE obj.ring256=",str(fpix[0])," or obj.ring256=",str(fpix[1])," or obj.ring256=",str(fpix[2])," or obj.ring256=",str(fpix[3])]) 
                        try:
                            dd=qc.query(sql=qtext,fmt="csv")
                        except Exception as e:
                            if "Query timeout" in e.message: #if the query timed out,
                                the_query=Table(names=qnames,dtype=qdtypes)
                                for iq in [0,1,2,3]:
                                    qqtext="".join(["SELECT meas.measid,meas.objectid,meas.exposure FROM nsc_dr1.meas as meas JOIN nsc_dr1.object as obj on obj.id=meas.objectid WHERE obj.ring256=",str(fpix[iq])]) #split query into 4, try again
                                    ddq=qc.query(sql=qtext,fmt="csv")
                                    the_query_temp=query_fix(ddq,qdtypes)
                                    the_query=vstack([the_query,the_query_temp])
                            else: print("query failed")
                        else: the_query=query_fix(dd,qdtypes)
                    # add the queried columns to the tracklet table!
                        t_out=join(t_out,the_query,keys="measid",join_type="left")

                # add exposure info to table
                    t_out=join(t_out,file_exposures,keys="exposure",join_type="left")
                # write a new fits file with just the tracklet mmts
                    #t_out.write("concat_files/hgroup_%s/cf_tracklets_hp_%s.fits" % (subdir,pix),format="fits",overwrite=True) #hyalite_msu
                    t_out.write(basedir+"hpix/hgroup_"+subdir+"/") #tempest_msu
                # add HP tracklet mmts to hgroup table
                    cat_mmts=vstack([cat_mmts,t_out])

                    # --Calculate info & write 80-col line for each tracklet observation--
                    for trl in np.unique(t_out['cluster_label']): #for every tracklet,
                        tr=t_out['cluster_label']==trl
                        tracklet=t_out[tr]
                        obs_num=0 #counter for tracklet observation number
                        for obs in tracklet: #for every observation,
                            # >calculate the necessary info
                            obs_80s=calc_mpc80(obs,scope_names,scope_codes) #obs_80s = [obs_date,obs_ra,obs_dec,obs_mag,obs_code]
                            # >create the line in 80-column format
                            linne = make_80col(obs_num,t_num,obs_80s[0],obs_80s[1],obs_80s[2],obs_80s[3],obs_80s[4])
                            # >write the line to hgroup MPC80 txt file
                            file80.write(linne)                             
                            obs_num+=1
                        t_num+=1
    file80.close()
    cat_mmts.write("concat_files/cat_files/canfind_mmts_hgroup_%s.fits"%subdir,overwrite=True)
    os.replace("canfind_tracklet_obs_%s.txt" % subdir, "concat_files/txt_files/canfind_tracklet_obs_%s.txt" % subdir)