#!/usr/bin/env python

# AUTHOR: Katie Fasbender
#         katiefasbender@montana.edu

# job_creator_v2.py is a python script that can submit jobs to the slurm job queue to your discretion.    
# Use if you have many jobs that you don't want to put in the slurm job queue all at once (this will 
# overload the queue - bad!).  This file was created to run CANFind on MSU's Tempest Research Cluster.  
# CANFind analyzes NSC measurements from 1 HEALPix (HP, NSIDE=128) at a time, source code in "canfind_v2.py".  
# This script should be run from the Tempest directory /home/x25h971/canfind_dr2.  

# CANFIND COMMAND: $python /path/to/canfind_v2.py <HP number> <analysis marker>
# This will analyze 1 healpix, with the assigned number <HPix number>, with CANFind.  

# This script writes a "job_name.sh" file for each job that will analyze a certain number "w" (#HP/job) 
# of a total "x" HP.  It will mantain "y" number of running jobs, checking every "z" seconds.

# "w" will depend on the number of measurements per HP (NMEAS), varying as we go through the HP list.   
# "x" depends on your HP input list; "x" and "y" are defined by you in the cmd line.  "z" is defined in the code.  

# COMMAND: $python /path/to/job_creator_v2.py /path/to/HP_list_.fits
# where job_creator_v2.py = this file (if you couldn't tell)
# and HP_list.fits = a fits file, please, with a list of HP to analyze.


#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------
import numpy as np
from argparse import ArgumentParser
from astropy.table import Table,Column
import os
import subprocess
import sys
import time 

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
def makedir(dir):
    '''makes a directory with name "dir" if it does not exist
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


def write_jscript(job_name,cmds):
    '''writes a SLURM job script to "job_name.sh"
       Lines starting with #SBATCH are read by Slurm. 
       Lines starting with ## are comments.
       All other lines are read by the shell
    Arguments:
    ----------
    job_name (str)
            name of job, job script file
    cmds (str list)
            python commands to run CF on HPs
    Returns:
    --------
    job_file (str)
            job filename the job script is written to
    '''
       	job_file = job_name+".sh"
        with open(job_file,'w') as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%d\n"%job_name)        # job name
            fh.writelines("#SBATCH --output=%d.out/n"%job_name)      # output file (%j = jobid)
            fh.writelines("#SBATCH --error=%d.err\n"%job_name)       # error file
            fh.writelines("#SBATCH --partition=unsafe\n")            # queue partition to run the job in
            fh.writelines("#SBATCH --ntasks=1\n")                    # for running in parallel
            fh.writelines("#SBATCH --nodes=1\n")                     # number of nodes to allocate
            fh.writelines("#SBATCH --ntasks-per-node 1\n")           # number of cores to allocate; set with care
            fh.writelines("#SBATCH --mem=60000\n")                   # memory, set --mem with care!!!!! refer to hyalite quickstart guide
            fh.writelines("#SBATCH --time=00:60:00\n")               # Maximum job run time
            fh.writelines("module load Anaconda3/2020.07\n")	     # load anaconda
            fh.writelines("source activate $HOME/condaenv/\n")       # load conda environment
            for cmd in cmds:
                fh.writelines(cmd.strip()+'\n')                      # write python CF command
        return job_file

#-----------------------------------------------------------------------------
# Main Code
#-----------------------------------------------------------------------------

if __name__ == "__main__":

    time0 = time.time() # start time

    # Initiate input arguments 
    #-------------------------
    parser = ArgumentParser(description='Run CANFind on selected HPix (NSIDE=128) from NSC DR2')
    parser.add_argument('--nrunjobs', type=int, nargs=1, default=1, help='Number of running jobs to maintain at any given time')
    parser.add_argument('--hpfile',type=str,nargs=1,default="healpix_good_dr2.fits",help='Filename of HEALPix list to analyze')
    args = parser.parse_args()

    # define inputs & variable values   
    nrunjobs = int(args.nrunjobs) # "y" number of jobs to maintian 
    sleep_time = 60               # "z" # of seconds to sleep between checking job channels
    hp_per_job = 5                # max number of HPix to analyze per job
    nmeas_tot = 10000             # max number of measurements per job, change definitely to smthn higher! 
    basedir = "/home/x25h971/canfind_dr2"
    hp_file = args.hplist         # full HP table  (input_list -> hp_file)
    if not os.path.exists(hp_file):                    
        print("%s doesn't exist!" %hp_file)
        sys.exit(1)
    hp_list = Table.read(hp_file)[0:20]                  #(hp_data -> hp_list) ***********REMOVE[0:20]***********
    nhp = len(hp_list) # total number of HP
    hp_list['cmd'] = Column(dtype="U100",length=nhp)
    hp_list['done'] = Column(np.zeros(length=nhp))
    hp_list['outfile'] = Column(dtype="U100",length=nhp)


    # Check HP list for completeness
    # ------------------------------
    for hpx in range(0,nhp_list):
        outdir = basedir+"/hpix/hgroup_"str(int(hp_list[hpx]['PIX'])//1000)
        makedir(outdir)
        outfile = outdir+"/healpix_"+str(hp_list[hpx]['PIX'])".fits"
        hp_list['outfile'][hpx] = outfile
        if not os.path.exists(outfile):
            hp_list['done'][hpx] = 1
        else:
            hp_list['cmd'][hpx] = "python "+basedir+"/files/canfind_v2.py "+str(hpx['PIX'])+" "+str(hpx['marker'])
    hp_torun = hp_list[hp_list['done']==0] # table of HP without outfiles (still must be analyzed) as of right now
    nhp_torun = len(hp_torun)              # total number "x" of HP to run 


    # Create job structure
    #---------------------
    jstruct_file = basedir+"/files/cf_jstruct_"+time0+".fits"
    jstruct = hp_torun.copy()
    jstruct['cputime'] = Column(dtype-"U100",length=nhp_torun)
    jstruct['jobchannel'] = Column(np.zeros(length=nhp_torun))
    jstruct['jobid'] = Column(np.repeat("-99.99",length=nhp_torun))
    jstruct['jobname'] = Column(dtype-"U100",length=nhp_torun)
    jstruct['jobnum'] = Column(np.zeros(length=nhp_torun))
    jstruct['jobstatus'] = Column(dtype="U100",length=nhp_torun)  
    jstruct['maxrss'] = Column(dtype-"U100",length=nhp_torun)
    jstruct['maxvmsize'] = Column(dtype-"U100",length=nhp_torun)
    jstruct['submitted'] = Column(np.zeros(length=nhp_torun))

    # loop through HP & parcel out among jobs 
    jstruct.sort('NMEAS')
    hpn = 0 # keep track of how many HP have been assigned
    jn = 0  # keep track of how many jobs have been "created"
    while hpn<=nhp_torun:
        # keep adding HP to a job until either:
        #(a) 5 HP in the job, or 
        #(b) NMEAS > nmeas_tot
        the_hps = [] # indices of HPs in job
        nm_tot = 0   # to keep track of #measurements in job 
        jobflag = 0
        while jobflag==0:
            the_hps.append(hpn)
            nmeas+=int(jstruct[hpn]['NMEAS'])
            if (len(the_hps)>3) or (nm_tot>nmeas_tot): 
                jstruct['jobnum'][the_hps] = jn
                jn+=1
                jobflag=1
            hpn+=1


    # Start submitting jobs 
    # ---------------------
    # loop through jobs 
    jb = 0  
    njobs=len(np.unique(jstruct['jobname'])) # total number of jobs to run
    while jb<=njobs:
        # loop through job channels
        for jbchannel in nrunjobs:

            # Check status of previous submitted job
            channel_ind = set(np.where(jstruct['jbchannel']==jbchannel)[0])
            submitted_ind = set(np.where(jstruct['submitted']==1)[0])
            unsubmitted_ind = set(np.where(jstruct['submitted']==0)[0])

            # get index & status of last job submitted 
            last_sub = list(channel_ind_ind & submitted_ind)
            if len(last_sub)==0: # if no jobs have been submitted to this channel
                lastjob=sort(list(channel_ind))[0]
            else: # else, grab the last job submitted
                lastjob=sort(last_sub)[-1] 
            last_jid = jstruct['jobid'][lastjob][0]
            if last_jid!="-99.99": lj_status = (subprocess.getoutput("sacct -n -X --format state --jobs="+last_jid).split("\n")[-1]).strip()
            else: lj_status = "NONE" #no jobs have been submitted
            jstruct['jobstatus'][lastjob] = lj_status
        
            # --If that job is still running or requeued or pending, pass ahead to sleep 
            if lj_status="PENDING" ls_status="REQUEUED" or lj_status="RUNNING":
                print("job "+str(last_jid)+" still running on "+str(jbchannel))

            # --Else, update statuses and write/submit a new job script
            else:
                # if last job was completed, get some info about it 
                if lj_status=="COMPLETED":
                    ljinfo = subprocess.getoutput("sacct -n -P --delimiter=',' --format cputimeraw,maxrss,maxvmsize --jobs "+last_jid)
                    ljinfo = ljinfo.split("\n")[-1].split(",")
                    jstruct['cputime'][lastjob] = ljinfo[0]
                    jstruct['maxrss'][lastjob] = ljinfo[1]
                    jstruct['maxvmsize'][lastjob] = ljinfo[2]
                # get indices of new job
                thisjob = (np.where(jstruct['jobnum']==jb)[0])
                # --Write & job script--
                job_name = "cf_dr2_"+jb          
                cmds = list(jstruct['cmd'][thisjob])
                job_file = write_jscript(job_name,cmds,basedir)
                # --Submit job script to slurm queue--
                os.system("sbatch %s" %job_file)
                jstruct['submitted'][thisjob] = 1
                # get jobid of submitted job
                jinfo = subprocess.getoutput("sacct -n -P --delimiter=',' --format jobid --name "+job_name)
                jinfo = jinfo.split("\n")[-1].split(",")
                jstruct['jobchannel'][thisjob] = jbchannel  
                jstruct['jobid'][thisjob] = jid
                jstruct['jobname'][thisjob] = job_name
                jb+=1

        # Save job structure, sleep before checking again   
        jstruct.write(jstruct_file,overwrite=True)
        time.sleep(sleep_time)  # SLEEP before checking/submitting next jobs 