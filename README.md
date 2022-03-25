# CANFind

CANFind (a Computationally Automated NSC-tracklet Finder) (sorry) was developed to detect "tracklets" (3+ measurements of a moving object) in the NOIRLab Source Catalog (NSC).  This repository contains the source code for CANFind v1 and v2, the code used to generate figures in the introductory article https://iopscience.iop.org/article/10.3847/1538-3881/ac2230 

CANFind was developed on the Astro Data Lab Jupyter notebook server.  It is designed to analyze one HEALPix (NSIDE=128,256) at a time, and takes approximately 5 seconds to 5 minutes for most HEALPix of typical object density.  Analysis is performed on Montana State University's Hyalite Computing Cluster (v1, NSC DR1) and Tempest Research Cluster (v2, NSC DR2), which allows for the analysis of multiple HP at a time, limited by the querying capacity of the NSC measurements & objects table. The python scripts used to generate job files for Hyalite's queue, managed by slurm, are also included in this repository. 

** A complete catalog of CANFind tracklets is available in FITS format at https://katiefasbender.wixsite.com/physics/research.


To run CANFind, you'll need...
- canfind.py
- NOIRLab's datalab query client, found here: https://github.com/noaodatalab-user/datalab-client 
- Healpy python package https://healpy.readthedocs.io/en/latest/install.html 
- scikit-learn python package https://scikit-learn.org/stable/install.html 

the "CANFind command" is:

```$ python path/to/canfind.py <HP_#> <analysis_marker> ```

where <analysis marker> is

  
To run CANFind on Hyalite or Tempest (for MSU members), you'll need...
- a list of HEALPix which you'd like to analyze, in the form of a FITS file (ne is provided here, "nsc_dr1_hp_canfind.fits") with 2 columns: 
    - "PIX" (HPix number, NSIDE=128)
    - "MARKER" (if 1: query the HPix, analyze mmts with CANFind, & save mmts and their tracklet labels; if 0: query the HPix & save the mmts not associated with SOs)
- job_creator.py 
    - you'll have to set some variables inside this file, guidelines in the comments 

the command is:

``` $ python path/to/job_creator.py path/to/<healpix_filename>.fits ```
