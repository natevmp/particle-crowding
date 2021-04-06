#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=24:0:0
#$ -wd ~/BrownianCrowding
#$ -j y
#$ -N BrownianParticlesManySims
#$ -o ~/BrownianCrowding
#$ -t 1-25

module load julia

DATADIR=$HOME/BrownianCrowding

rsync -rltv $DATADIR/ $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

# julia HPC/addPackages.jl
julia HPC/growingPop_script.jl ${SGE_TASK_ID} 0.005 4500

# rsync -rltv $TMPDIR/ $DATADIR/
mv growingPop_multSims*.jld $DATADIR/

echo "success"

