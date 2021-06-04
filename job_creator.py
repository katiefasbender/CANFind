#!/usr/bin/env python

# AUTHOR: Katie Fasbender
# 	  katiefasbender@montana.edu

# job_creator.py is a python script that can submit "batches" of jobs to slurm (on Hyalite) at your discretion!  
# Use if you have many jobs that you don't want to put in the slurm job queue all at once (this will overload the queue - bad!)

# This file was created to run CANFind on Hyalite.  CANFind analyzes NSC measurements from 1 HEALPix (HP, nside=128) at a time.
# The "CANFind command": python canfind.py <HPix number> <analysis marker>
# The CANFind command will send 1 healpix, with the assigned number <HPix number>, through CANFind.  

# Each CANFind job will analyze a set number "x" of HP.
# For each job, this script writes a "job_name.sh" file.
# In each "job_name.sh" file created by this script, "x" lines of CANFind commands are written. 
# "y" number of job files will be submitted in a batch before sleeping for "z" seconds.

# input formula (write in command line to run this file): python path/to/job_creator.py path/to/HP_list_filename.fits
# where job_creator.py = this file (if u couldn't tell) which will write and submit a batch "y" number of jobs every "z" seconds to slurm
# and HP_list_filename.fits = a fits file with a list of healpix to analyze.


#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
from astropy.io import fits
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

	input_file=sys.argv[1] #grabs "path/to/HP_list_filename.fits"
	if not os.path.exists(input_file): #to check whether there actually is a file.
		print("You messed up!!!!!!!  %s doesn't exist!" %input_file)
		sys.exit(1)

	hdul=fits.open(input_file)
	data=hdul[1] #the list of HP that must be analyzed, from the fits file

## the following is in the "job_name.sh" file.
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.

	#set some parameters:
	hp_per_job=1 #"x" number of HP to be analyzed per job
	num_hp=len(hdul[1].data) #total # of HP to be analyzed
	max_range=int(num_hp/hp_per_job) #total number of jobs
	max_jobs=1 #"y" number of jobs in a batch to be submitted to slurm queue before "sleeping" (this avoids overloading slurm)
	sleep_time=300 #"z" number of seconds to sleep between job batches
	m=0 #a counter to keep track of how many job files have been written in current batch

	for jb in range(0,max_range): #for each job, let's write a job file
		job_file='job_name_%d.sh' % jb #define the name of the job file about to be written ("job_name.sh")
		with open(job_file,'w') as fh: #even if the file doesn't exist, it will be created and all these lines written to it
			fh.writelines("#!/bin/bash\n")
			fh.writelines("#SBATCH --job-name=name_of_job_%j")       # job name
			fh.writelines("#SBATCH --output=nft-%j.out\n")           # standard output file (%j = jobid)
			fh.writelines("#SBATCH --error=nft-%j.err\n")            # standard error file
			fh.writelines("#SBATCH --partition=unsafe\n")	         # queue partition to run the job in
			fh.writelines("#SBATCH --ntasks=1")                      #for running in parallel
			fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
			fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care 
			fh.writelines("#SBATCH --mem=200\n")                     # memory, set --mem with care!!!!! refer to hyalite quickstart guide
			fh.writelines("#SBATCH --time=00:60:00\n")               # Maximum job run time
			fh.writelines("module load Anaconda3/5.1.0\n")           #load anaconda, needed for running python on Hyalite!
			for i in range(0,hp_per_job):#for every HP included in this job, 
				counter=i+(hp_per_job*jb) #to keep track of how many HP have been submitted
				num_pix=hdul[1].data[counter][0] #the HP_#
				pix_marker=hdul[1].data[counter][1] #the marker (1 for CANFind, 0 for...NOT CANFind)
				subdir=int(int(num_pix)//1000)
				makedir("canfind_hpix/hgroup_%d" % subdir)
				if counter<(num_hp+1): #if we haven't exceeded the total number of HP to analyze,
					fh.writelines("python canfind.py %d %d\n" % (num_pix,pix_marker)) #write the CANFind command to the job file 		
		os.system("sbatch %s" %job_file) #send the job file to slurm

		#check to see if the batch of files has been submitted to slurm
		# if so, time to sleep so we don't overload the slurm queue.  
		m=m+1  
		if m==max_jobs: #after "max_jobs" number of jobs (set above), sleep for "sleep_time" number of seconds before submitting more jobs
			time.sleep(sleep_time) #SLEEP COMMAND   
			m=0 #reset the counter 
