#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=48:0:0
#$ -wd ~/BrownianCrowding
#$ -j y
#$ -N BrownianParticlesManySims
#$ -o ~/BrownianCrowding
#$ -t 1-6

module load julia

DATADIR=$HOME/BrownianCrowding

rsync -rltv $DATADIR/ $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

# julia HPC/addPackages.jl
julia HPC/fixedDensity_script.jl ${SGE_TASK_ID}

# rsync -rltv $TMPDIR/ $DATADIR/
mv simResultFixedDensity_*.jld2 $DATADIR/

echo "success"

