#!/bin/bash 

#SBATCH -t 8640 # Runtime in minutes 
#SBATCH --mem=32000 # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH --mem-per-cpu=32000 # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH --open-mode=append 
#SBATCH --output=32core_%j.out # Standard out goes to this file
#SBATCH --nodes=1 # Ensure that all cores are on one machine (for IDL license reasons)
#SBATCH --ntasks=32 # request ntask/ntasks-per-core cores
#SBATCH --ntasks-per-core=1

# load IDL
module load IDL/8.5.1-fasrc01

# get environment set up
source /n/home12/jeastman/.exofast.tcshrc

# one license used per node, only 62 licenses available
# but we can run one fit per core (32 cores per node) 
for i in `seq 32`; do idl -e 'fittess' & sleep 60; done

# wait for all IDL tasks to finish
for pid in `ps -ef | grep idl | awk '{print $2}' | head -n -1`; do wait ${pid}; done

## run like this (use 50/62 licenses)
# for i in `seq 50`; do sbatch fit.sbatch; done

### useful commands ###
## show the scheduling queue 
# showq-slurm -p serial_requeue -o

## show the status of job 34002694
#squeue -j 34002694

## display the number of IDL licenses in use
# /n/sw/idl-8.5.1/idl85/bin/lmstat -a
