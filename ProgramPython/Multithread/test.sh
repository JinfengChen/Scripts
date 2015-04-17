#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

python test.py

echo "Done"
