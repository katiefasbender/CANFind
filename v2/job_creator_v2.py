#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# job_creator_v2.py is a python script that can submit "batches" of jobs to slurm (on Hyalite) at your discretion!  
# Use if you have many jobs that you don't want to put in the slurm job queue all at once (this will overload the queue - bad!)
# This file was created to run CANFind on Hyalite.  CANFind analyzes NSC measurements from 1 HEALPix (HP, NSIDE=128) at a time.
# It is meant to be run from the directory /scratch/katie.fasbender/canfind_dr2_w and
# outfiles are written to /scratch/katie.fasbender/canfind_hpix_dr2

# CANFIND COMMAND: python canfind_v2.py <HPix number> <analysis marker>
# This will analyze 1 healpix, with the assigned number <HPix number>, with CANFind.  

# This script writes a "job_name.sh" file for each job that will analyze "x" number of HP.
# In each "job_name.sh" file, "x" lines of CANFind commands are written. 
# "y" number of job files will be submitted in a batch before sleeping for "z" seconds.

# INPUT FORMULA (write in command line to run this file): python <path/to/job_creator_v2.py> <path/to/HP_list_filename.fits>
# where job_creator_v2.py = this file (if u couldn't tell) which will write and submit a batch of "y" number of jobs every "z" seconds to slurm
# and HP_list_filename.fits = a fits file with a list of healpix to analyze.


#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
from astropy.table import Table,Column
import sys
import os
import time 

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def makedir(dir):
    '''makes a directory with name "dir" if it does not exist'''
    if not os.path.exists(dir):
        os.mkdir(dir)

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    # Get the full list of HP 
    input_file=sys.argv[1]                                # grabs "path/to/HP_list_filename.fits"
    if not os.path.exists(input_file):                    # to check whether there actually is a file.
        print("%s doesn't exist!" %input_file)
        sys.exit(1)
    hp_data=Table.read(input_file)                        # table with full HP list 

    # Get the list of HP with columns for status & submission (1/0)
    hp_done_file="done_"+input_file
    if not os.path.exists(hp_done_file):                 
        done_data=hp_data.copy()
        done_data.add_column(Column(np.zeros(len(hp_data)),name="status"))
        done_data.add_column(Column(np.repeat(-1,len(hp_data)),name="job_num"))
    else: done_data=Table.read(hp_done_file)              # table with HP outfile statuses
    
    data=done_data[done_data['status']==0]                # table of HP without outfiles (still must be analyzed) as of right now

#-------------------
# Writing job script
#--------------------------------------------------------------------------------------------------------------------------------- 
# The following code writes lines to the "job_name.sh" file.
# Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
# All other lines are read by the shell.

    # Set some parameters:
    hp_per_job=2#1                                  # "x" number of HP to be analyzed per job
    num_hp=len(data)                              # total # of HP to be analyzed
    max_range=int(num_hp/hp_per_job)              # total number of jobs
    max_jobs=4                                    # "y" # of jobs in a batch to be submitted to slurm queue before "sleeping" (to avoid overloading slurm)
    sleep_time=60#300                                # "z" # of seconds to sleep between job batches
    y=0                                           # counter to keep track of (# job files) written in current batch
    outfile_array=[[] for _ in range(0,max_jobs)]    # an array to store names of all outfiles that should be created by this job batch
    increm=0
    # For each job "jb", write a "job_name.sh" file
    for jb in range(0,max_range):
        job_file='cfdr2_%d.sh' % jb                    # define the name of the job file about to be written ("job_<jb>.sh")
        with open(job_file,'w') as fh:               # even if the file doesn't exist, it will be created and all these lines written to it
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=cfdr2_%d\n"%jb)     # job name
            fh.writelines("#SBATCH --output=cfdr2_%d"%(jb)+"_%j.out\n")           # standard output file (%j = jobid)
            fh.writelines("#SBATCH --error=cfdr2_%d.err\n"%jb)            # standard error file
            fh.writelines("#SBATCH --partition=unsafe\n")	     # queue partition to run the job in
            fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
            fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
            fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care 
            fh.writelines("#SBATCH --mem=2000\n")                     # memory, set --mem with care!!!!! refer to hyalite quickstart guide
            fh.writelines("#SBATCH --time=00:60:00\n")               # Maximum job run time
            fh.writelines("module load Anaconda3/5.1.0\n")           # load anaconda, needed for running python on Hyalite!

            # Loop through appropriate number of hpix per job
            i=0
            counter=(hp_per_job*jb)+increm                                         # to keep track of how many HP have been submitted
            while i<hp_per_job:
            #for i in range(0,hp_per_job):
                num_pix=data[counter]['PIX']                                      # <HP_#>
                pix_marker=data[counter]['marker']                                # the marker (1 for CANFind, 0 for...NOT CANFind)
                subdir=int(int(num_pix)//1000)                                    # the subdirectory number to write outfile to 
                outdir="/mnt/lustrefs/scratch/katie.fasbender/canfind_hpix_dr2/hgroup_"  
                makedir(outdir+str(subdir))
                outfile_name="healpix_%s.fits" % str(num_pix)

                # Write the CANFind command to jobfile, if (a) there are still hpix to be analyzed 
                #                                          (b the HP outfile does not already exist
                if (counter<(num_hp+1)):                   # check for case (a)
                    hpind=(done_data['PIX']==num_pix).tolist().index(True)
                    if not os.path.exists(outdir+str(subdir)+"/"+outfile_name): #check for case(b)
                        print("hpix to be analyzed = ",num_pix)
                        fh.writelines("python canfind_v2.py %d %d\n" % (num_pix,pix_marker))
                        done_data[hpind]['job_num']=jb
                        outfile_array[y].append(outdir+str(subdir)+"/"+outfile_name)      
                        i+=1
                        counter+=1
                    #elif os.path.exists(outdir+str(subdir)+"/"+outfile_name):
                    else:
                        done_data[hpind]['status']=1
                        counter+=1
                        increm+=1
                else: i=hp_per_job
        done_data.write(hp_done_file,overwrite=True)
                        
#-------------------
# Job script written
#--------------------------------------------------------------------------------------------------------------------------------- 

        # Submit the jobfile to slurm queue
        os.system("sbatch %s" %job_file)
#---------------------
# Job script submitted
#--------------------------------------------------------------------------------------------------------------------------------- 

        # Check to see if the maximum # of jobfiles has been submitted to slurm.  
        # If so, time to sleep before submitting the next batch of jobfiles 
        # so we don't overload the slurm queue and the datalab querying capabilities.

        y=y+1                              # increment the counter to indicate submitting last job
        if y==max_jobs:                    # once the maximum # of jobs are currently submitted & running,
            while y==max_jobs:             # every sleep_time seconds, check which outfiles have been written
                print("sleeping for ",sleep_time," seconds")
                time.sleep(sleep_time)                            # SLEEP
                print("checking for outfiles...")
                job_num_outfiles=[]                               # an array for the number of files left (each element for 1 job)
                
                for jnum in range(0,max_jobs):                    # for every job in the batch, check the outfiles
                    outfiles_array_cp=outfile_array[jnum].copy()   # the latest list of outfiles to check for this job
                    hpjid=[]
                    #print("outfiles of this job = ",outfiles_array_cp)
                    for ofile in outfiles_array_cp:               # for each outfile, check if exists.  If so, remove from array
                        ofile_hpix=int((ofile.split("_")[-1]).split(".")[0])
                        hpjid.append(done_data[done_data['PIX']==ofile_hpix]['job_num'])
                        if os.path.exists(ofile): 
                            outfile_array[jnum].remove(ofile)
                            hpixid=(done_data['PIX']==ofile_hpix).tolist().index(True)
                            done_data[hpixid]['status']=1
                    #print("hpjid = ",hpjid)
                    with open("cfdr2_%d.err"%hpjid[0]) as jf:      # check the error file for the job status
                        jflines=jf.readlines()
                    if len(jflines)!=0:                             # if there was an error, call it a failed job.
                        outfile_array[jnum]=[]
               
                    job_num_outfiles.append(len(outfile_array[jnum]))
                
                done_data.write(hp_done_file,overwrite=True)
                jno=np.array(job_num_outfiles)
                job_zeros=jno[jno==0]
                outfile_array=[_ for _ in outfile_array if _!=[]]
                outfile_array=outfile_array+[[] for _ in range(0,len(job_zeros))]
                print("files/job = ",job_num_outfiles)
                print(outfile_array)
                #print("y = ",y-len(job_zeros))
                y-=int(len(job_zeros))

#-----------------------------
# "y" Jobs maintained in batch 
#--------------------------------------------------------------------------------------------------------------------------------- 
