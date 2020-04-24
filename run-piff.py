#!/usr/bin/env python
"""
Piff job submission script
"""
__author__ = "Alex Drlica-Wagner"
import os
import numpy as np
import subprocess
import glob

import argparse
parser = argparse.ArgumentParser(description=__doc__)
group = parser.add_mutually_exclusive_group()
group.add_argument('--sbatch',action='store_true',
                   help="submit with sbatch from login node")
group.add_argument('--srun',action='store_true',
                   help="submit with srun from interactive node")
parser.add_argument('-d','--dryrun',action='store_true',
                    help="generate submission files but don't submit")
parser.add_argument('-f','--force',action='store_true',
                    help="force rerun")
parser.add_argument('-r','--runid',default=0,type=int,
                    help="integer run id for output directory")
parser.add_argument('-b','--bands',action='append',nargs='+',
                    help="select specific bands")
group = parser.add_mutually_exclusive_group()
group.add_argument('-v','--visit',default=None,action='append',
                    help="specify visit number(s)")
group.add_argument('-n','--nvisit',default=None,type=int,
                    help="number of visits to run")
parser.add_argument('-j','--jobs_per_slurm',default=10,type=int,
                    help="number of visits to run per slurm job")
parser.add_argument('--vmin',default=None,type=int,
                    help="minimum visit number")
parser.add_argument('--vmax',default=None,type=int,
                    help="maximum visit number")
args = parser.parse_args()

#args.dryrun=False # True #hack
runid = args.runid

basedir  = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1'
outbase   = '/global/cscratch1/sd/mjarvis/DC2/dc2.1'
config   = os.path.join(outbase,'piff.yaml')

dirnames = sorted(glob.glob(os.path.join(basedir,'calexp/*')))
visitid  = np.char.rpartition(dirnames,'/')[:,-1]
visit   = np.char.rpartition(visitid,'-')[:,0].astype(int)
band    = np.char.rpartition(visitid,'-')[:,-1]

alldata = np.rec.fromarrays([visitid,visit,band],names=['visitid','visit','band'])
sel = np.ones(len(alldata),dtype=bool)

if args.bands:
    sel &= np.in1d(alldata['band'],args.bands)

if args.visit:
    sel &= np.in1d(alldata['visit'],args.visit)
else:
    if args.vmin:
        sel &= alldata['visit'] >= args.vmin
    if args.vmax:
        sel &= alldata['visit'] <= args.vmax

data = alldata[sel]

# old
#SBATCH -C knl
#SBATCH --qos debug
#SBATCH --cpus-per-task=16

# Template for the piff command
piff = "piffify -l {logfile} {config} output.dir={outdir} input.image_file_name={image_file_name} output.file_name={file_name} output.stats.0.file_name={rho_name} output.stats.1.file_name={shape_name} output.stats.2.file_name={twod_name} output.stats.3.file_name={outcat_name}"

# Template for slurm
slurm = """#!/bin/bash -l
#SBATCH --account=m1727
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=04:00:00
#SBATCH --job-name={visit}
#SBATCH --output={logfile}

export OMP_NUM_THREADS=8
source /usr/common/software/python/2.7-anaconda-2019.07/etc/profile.d/conda.sh
conda activate piff
"""

run = "{cmd} >& {logfile}"
srun = "srun -N 1 {cmd} >& {logfile} &"

visit_list = []
complete = 0

# First figure out which visits we're going to run.
for i,d in enumerate(data):
    if args.nvisit and len(visit_list) >= args.nvisit:
        break

    visit = d['visitid']
    if args.visit and not np.in1d(visit,args.visit):
        print("Skipping {visit}...".format(visit=visit))
        continue

    outdir=os.path.join(outbase,'piff-run-{runid:04d}/{visit}'.format(runid=runid,visit=visit))
    outcat_name = 'piff_{visit}_cat.fits'.format(visit=visit)
    outfile = os.path.join(outdir,outcat_name)
    if os.path.exists(outfile) and not args.force:
        print("Found {}; skipping...".format(outfile))
        complete += 1
        continue

    visit_list.append(visit)
print('{} of {} visits are already complete.'.format(complete,len(data)))
print('Starting jobs to run {} visits.'.format(len(visit_list)))

# Now run them.
for i,visit in enumerate(visit_list):
    print("Visit",visit)

    outdir=os.path.join(outbase,'piff-run-{runid:04d}/{visit}'.format(runid=runid,visit=visit))
    file_name = 'piff_{visit}.fits'.format(visit=visit)
    outfile = os.path.join(outdir,file_name)

    image_file_name="calexp/{visit}/*/calexp_{visit}-*.fits".format(visit=visit)
    run_outdir=os.path.join(outbase,'piff-run-{runid:04d}'.format(runid=runid))
    rho_name = 'piff_{visit}_rho.png'.format(visit=visit)
    shape_name = 'piff_{visit}_shape.png'.format(visit=visit)
    twod_name = 'piff_{visit}_twod.png'.format(visit=visit)
    outcat_name = 'piff_{visit}_cat.fits'.format(visit=visit)
    logfile = 'piff_{visit}.log'.format(visit=visit)

    if i % args.jobs_per_slurm == 0:
        slurm_logfile = os.path.join(run_outdir,logfile)
        subfile = os.path.join(run_outdir,'piff_{visit}.sub'.format(visit=visit))
        with open(subfile,'w') as sub:
            sub.write(slurm.format(visit=visit,logfile=slurm_logfile))

    # Make the output directory
    subprocess.call('mkdir -p {outdir}'.format(outdir=outdir),shell=True)

    # Create piff command
    logfile = os.path.join(outdir,logfile)
    cmd = piff.format(config=config, outdir=outdir, image_file_name=image_file_name,
                      file_name=file_name, rho_name=rho_name, shape_name=shape_name,
                      twod_name=twod_name, outcat_name=outcat_name, logfile=logfile)
    env = dict(os.environ, OMP_NUM_THREADS='2')

    if args.sbatch:
        # Write the submission script
        submit="sbatch {subfile}".format(subfile=subfile)
        exe = submit
        with open(subfile,'a') as sub:
            sub.write("srun {cmd}\n".format(cmd=cmd))
        # Run this if the next job will start a new submission script
        # Or if we are on the last visit.
        run_job = (i+1) % args.jobs_per_slurm == 0 or i == len(visit_list)-1
    elif args.srun:
        exe = srun.format(cmd=cmd,logfile=logfile)
        run_job = True
    else:
        exe = run.format(cmd=cmd,logfile=logfile)
        run_job = True

    print(exe)

    if not args.dryrun and run_job:
        print("Executing job...")
        subprocess.call(exe, env=env, shell=True)
