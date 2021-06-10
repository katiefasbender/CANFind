# CANFind

CANFind (a Computationally Automated NSC-tracklet Finder) (sorry) was developed to detect "tracklets" (3+ measurements of a moving object) in the NOIRLab Source Catalog.  This repository contains the source code for CANFind V1, the code used to generate figures in the paper submitted by CANFind's creators,  

CANFind was developed on the Astro Data Lab Jupyter notebook server.  It is designed to analyze HEALPix (NSIDE=128) at a time, and takes approximately 5 seconds to 5 minutes for most HEALPix of typical object density.  Analysis is performed on Montana State University's Hyalite Computing Cluster, which allows for the analysis of up to 10 HP at a time, limited by the querying capacity of the NSC measurements & objects table. The python scripts used to generate job files for Hyalite's queue, managed by slurm, are also included in this repository. 


To run CANFind, you'll need...
- canfind.py
- NOIRLab's datalab query client, found here: https://github.com/noaodatalab-user/datalab-client 

the "CANFind command" is:
$ python path/to/canfind.py <HPix #> <analysis marker>

  
To run CANFind on Hyalite (for MSU members), you'll need...
- a Hyalite account 
- a list of HEALPix which you'd like to analyze, in the form of a FITS file (ne is provided here, "nsc_dr1_hp_canfind.fits") with 2 columns: 
    - "PIX" (HPix number, NSIDE=128)
    - "MARKER" (if 1: query the HPix, analyze mmts with CANFind, & save mmts and their tracklet labels; if 0: query the HPix & save the mmts not associated with SOs)
- job_creator.py 
    - you'll have to set some variables inside this file, guidelines in the comments 

the command is:
$ python path/to/job_creator.py path/to/<healpix_filename>.fits
  
  
  
  
  
  
Tracklet Catalog Description:
  
  Column Names | Units | dtype | Description
               |       |       |
  'ra'         | deg   | float | mean tracklet RA
  'dec'        | deg   | float | mean tracklet dec
  'vra'        | ''/hr | float | tracklet proper motion in RA
  'vraerr'     | ''/hr | float | error in vra 
  'vdec'       | ''/hr | float | tracklet proper motion in Dec
  'vdecerr'    | ''/hr | float | error in vdec
  'el'         | deg   | float | mean tracklet ecliptic longitude l
  'eb'         | deg   | float | mean tracklet ecliptic latitude b
  'vel'        | ''/hr | float | tracklet proper motion in l
  'veb'        | ''/hr | float | tracklet proper motion in b
  'mjd'        | day   | float | first detection of tracklet in Modified Julien Date 
  'd_mjd'      | day   | float | time between first and last detection of tracklet
  'mag_auto'   | mag   | float | mean tracklet magnitude across all bands 
  'nmeas'      |       | int   | number of measurements in tracklet
  'tracklet_id'|       | int   | unique tracklet ID 
  'pix128'     |       | int   | tracklet HEALPix number (NSIDE=128)
  'redchi'     |       | float | reduced chi squared of tracklet on proper motion
  "u"          | mag   | float | mean tracklet magnitude in u band 
  "g"          | mag   | float | mean tracklet magnitude in g band
  "r"          | mag   | float | mean tracklet magnitude in r band
  "i"          | mag   | float | mean tracklet magnitude in i band
  "z"          | mag   | float | mean tracklet magnitude in z band
  "Y"          | mag   | float | mean tracklet magnitude in Y band
  "VR"         | mag   | float | mean tracklet magnitude in VR band 
  "pm_tot"     | ''/hr | float | total proper motion of tracklet
  "pmerr"      | ''/hr | float | uncertainty in proper motion

