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
group.add_argument('--sbatch',action='store_true')
group.add_argument('--srun',action='store_true')
parser.add_argument('-f','--force',action='store_true',
                    help="force rerun")
parser.add_argument('-r','--runid',default=0,type=int,
                    help="run id")
parser.add_argument('-v','--visit',default=None,action='append',
                    help="specify visit number(s)")
parser.add_argument('-b','--bands',action='append',nargs='+',
                    help="select specific bands")
parser.add_argument('-n','--nvisit',default=None,type=int,
                    help="number of visits to run")
parser.add_argument('--vmin',default=None,type=int,
                    help="minimum visit number")
parser.add_argument('--vmax',default=None,type=int,
                    help="maximum visit number")
parser.add_argument('-d','--dryrun',action='store_true',
                    help="generate submission files but don't submit")
args = parser.parse_args()

#args.dryrun=False # True #hack
runid = args.runid

basedir  = '/global/cscratch1/sd/desc/DC2/data/Run2.1i/rerun/calexp-v1'
outbase   = '/global/cscratch1/sd/mjarvis/DC2/dc2.1_test'
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
piff = "piffify {config} output.dir={outdir} input.image_file_name={image_file_name} output.file_name={file_name}"

# Template for slurm
slurm = """#!/bin/bash -l
#SBATCH --account=m1727
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --job-name={visit}
#SBATCH --output={logfile}

source /usr/common/software/python/2.7-anaconda-2019.07/etc/profile.d/conda.sh
conda activate piff
srun {cmd}
"""

run = "{cmd} >& {logfile}"
srun = "srun -N 1 {cmd} >& {logfile} &"

nsub = 0
for i,d in enumerate(data):
    if args.nvisit and nsub > args.nvisit:
        break

    visit = d['visitid']
    if args.visit and not np.in1d(visit,args.visit):
        print("Skipping {visit}...".format(visit=visit))
        continue

    image_file_name="calexp/{visit}/*/calexp_{visit}-*.fits".format(visit=visit)
    outdir=os.path.join(outbase,'piff-run-{runid:04d}/{visit}'.format(runid=runid,visit=visit))

    file_name = 'piff_{visit}.fits'.format(visit=visit)
    outfile = os.path.join(outdir,file_name)
    logfile = os.path.join(outdir,'piff_{visit}.log'.format(visit=visit))
    subfile = os.path.join(outdir,'piff_{visit}.sub'.format(visit=visit))

    if os.path.exists(outfile) and not args.force:
        print("Found {}; skipping...".format(outfile))
        continue

    # Make the output directory
    subprocess.call('mkdir -p {outdir}'.format(outdir=outdir),shell=True)

    # Create piff command
    cmd = piff.format(config=config,outdir=outdir,image_file_name=image_file_name,file_name=file_name)

    if args.sbatch:
        # Write the submission script
        submit="sbatch {subfile}".format(subfile=subfile)
        with open(subfile,'w') as sub:
            sub.write(slurm.format(visit=visit,logfile=logfile,cmd=cmd))
        exe = submit
    elif args.srun:
        exe = srun.format(cmd=cmd,logfile=logfile)
    else:
        exe = run.format(cmd=cmd,logfile=logfile)

    print(exe)

    if not args.dryrun:
        #subprocess.call(submit,shell=True)
        print("Executing job...")
        subprocess.call(exe,shell=True)

    nsub += 1
