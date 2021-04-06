#!/bin/bash -l
#PBS -l nodes=1:ppn=36,walltime=8:00:00
#PBS -A lp_nanomotors
#PBS -N many_gp

set -x

module load Julia

JOB_ID=$(echo $PBS_JOBID | awk -F. '{print $1}')

cd ${VSC_SCRATCH_NODE}

cp -r $HOME/particle-crowding ./
cd particle-crowding/test

parallel -j 36 julia growingPop_script.jl {1} ${RHO} ::: {a..b}{a..b}

mkdir ${VSC_SCRATCH}/many_gp_${JOB_ID}

mv growingPop_multSims_rho*.jld2 ${VSC_SCRATCH}/many_gp_${JOB_ID}

