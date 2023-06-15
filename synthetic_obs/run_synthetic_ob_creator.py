"""
Master Script that Submits Jobs to Create Synthetic Obs

shawn.s.murdzek@noaa.gov
Date Created: 15 June 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import os
import sys
import datetime as dt
import time
import pandas as pd


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# BUFR time stamps
bufr_time = [dt.datetime(2022, 4, 29, 12+i) for i in range(3)]

# Parameters passed to create_conv_obs.py for each prepBUFR time
nfiles = len(bufr_time)
wrf_dir = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/UPP'] * nfiles
bufr_dir = ['/work2/noaa/wrfruc/murdzek/real_obs/obs_rap_csv'] * nfiles
fake_bufr_dir = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/conv/data'] * nfiles
log_dir = ['/work2/noaa/wrfruc/murdzek/nature_run_spring/synthetic_obs_csv/conv/log'] * nfiles

# BUFR tags
bufr_tag = ['rap', 'rap_e', 'rap_p']

# Number of simultaneous jobs that can be running at once
max_jobs = 25
user = 'smurdzek'
allocation = 'wrfruc'
job_csv_name = 'synthetic_ob_creator_jobs.csv'


#---------------------------------------------------------------------------------------------------
# Create and Submit Job Scripts
#---------------------------------------------------------------------------------------------------

# First determine how many jobs are running and which jobs have already been submitted
njobs = len(os.popen('squeue -u %s' % user).read().split('\n')) - 2
if os.path.isfile(job_csv_name):
    sub_jobs = pd.read_csv(job_csv_name)
else:
    tmp_dict = {'date':[], 'idate':[], 'tag':[], 'submitted':[]}
    for i, t in enumerate(bufr_time):
        for tag in bufr_tag:
            tmp_dict['date'].append(t)
            tmp_dict['idate'].append(i)
            tmp_dict['tag'].append(tag)
            tmp_dict['submitted'].append(False)
    sub_jobs = pd.DataFrame(tmp_dict)

# Submit new jobs
while njobs < max_jobs:

    # End script if all jobs have been submitted
    if np.all(sub_jobs['submitted']):
        break

    idx = sub_jobs.loc[sub_jobs['submitted'] == False, 'idate'].values[0]
    tag = sub_jobs.loc[sub_jobs['submitted'] == False, 'tag'].values[0]

    # Determine first and last WRF file time
    wrf_start = None
    hr = 3
    while wrf_start == None:
        tmp = bufr_time[idx] - dt.timedelta(hours=hr)
        if os.path.isfile('%s/%s' % (wrf_dir[idx], tmp.strftime('%Y%m%d/wrfnat_%Y%m%d%H%M.grib2'))):
            wrf_start = tmp
            break
        hr = hr + 1

    tmp = bufr_time[idx] + dt.timedelta(hours=1)
    if os.path.isfile('%s/%s' % (wrf_dir[idx], tmp.strftime('%Y%m%d/wrfnat_%Y%m%d%H%M.grib2'))):
        wrf_end = tmp
    else:
        wrf_end = bufr_time[idx]

    # Create job script
    j_name = '%s/syn_obs_%s_%s.sh' % (log_dir[idx], bufr_time[idx].strftime('%Y%m%d%H%M'), tag)
    fptr = open(j_name) 
    fptr.write('#!/bin/sh\n\n')
    fptr.write('#SBATCH -A %s\n' % allocation)
    fptr.write('#SBATCH -t 08:00:00\n')
    fptr.write('#SBATCH --nodes=1 --ntasks=1\n')
    fptr.write('#SBATCH --mem=30GB\n')
    fptr.write('#SBATCH -o %s.%s.log\n' % (bufr_time[idx].strftime('%Y%m%d%H%M'), tag))
    fptr.write('#SBATCH --partition=orion\n\n')
    fptr.write('date\n. ~/.bashrc\nmy_py\n\n')
    fptr.write('python create_conv_obs.py %s \ \n' % wrf_dir[idx])
    fptr.write('                          %s \ \n' % bufr_dir[idx])
    fptr.write('                          %s \ \n' % fake_bufr_dir[idx])
    fptr.write('                          %s \ \n' % bufr_time[idx].strftime('%Y%m%d%H'))
    fptr.write('                          %s \ \n' % wrf_start.strftime('%Y%m%d%H'))
    fptr.write('                          %s \ \n' % wrf_end.strftime('%Y%m%d%H'))
    fptr.write('                          %s \ \n\n' % tag)
    fptr.write('date')
    fptr.close()
   
    # Submit job and mark it as submitted
    os.system('sbatch %s' % j_name)
    sub_jobs.loc[idx, 'submitted'] = True

    # Wait a few seconds for the job to submit, then check number of jobs
    time.sleep(15)
    njobs = len(os.popen('squeue -u %s' % user).read().split('\n')) - 2

# Save updated job submission CSV
sub_jobs.to_csv(job_csv_name)


"""
End run_synthetic_ob_creator.py
"""
