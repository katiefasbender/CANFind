# CANFind
CANFind

CANFind (a Computationally Automated NSC-tracklet Finder) (sorry) was developed to detect "tracklets" (3+ measurements of a moving object) in the NOIRLab Source Catalog.  This repository contains the source code for CANFind V1, the code used to generate figures in the paper submitted by CANFind's creators,  

CANFind was developed on the Astro Data Lab Jupyter notebook server.  It is designed to analyze HEALPix (NSIDE=128) at a time, and takes approximately 5 seconds to 5 minutes for most HEALPix of typical object density.  Analysis is performed on Montana State University's Hyalite Computing Cluster, which allows for the analysis of up to 10 HP at a time, limited by the querying capacity of the NSC measurements & objects table. The python scripts used to generate job files for Hyalite's queue, managed by slurm, are also included in this repository.  
