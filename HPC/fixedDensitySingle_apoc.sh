#!/bin/bash
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -l h_rt=48:0:0
#$ -wd ~/BrownianCrowding
#$ -j y
#$ -N BrownianParticlesManySims
#$ -o ~/BrownianCrowding

module load julia

DATADIR=$HOME/BrownianCrowding

rsync -rltv $DATADIR/ $TMPDIR/

cd $TMPDIR

echo "starting julia script..."

julia HPC/addPackages.jl
julia HPC/fixedDensity_script.jl 6

# rsync -rltv $TMPDIR/ $DATADIR/
mv simResultFixedDensity_*.jld2 $DATADIR/

echo "success"

