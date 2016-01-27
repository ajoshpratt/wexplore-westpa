#!/bin/bash -l
#$ -N CH4.Eq.14
#$ -pe fill_up 24
#$ -l h_rt=96:00:00
#$ -q compute48.q
#$ -V
#$ -cwd

source env.sh || exit 1

# start server

$WEST_ROOT/bin/w_run --work-manager=processes --n-workers=6 &> west.log &

wait
