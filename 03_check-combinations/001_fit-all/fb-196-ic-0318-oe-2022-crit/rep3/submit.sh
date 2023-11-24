#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J master_1_196
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  master_1_196-%J.o
#BSUB -e  master_1_196-%J.e
#
# Set any cpu options.
#BSUB -M 2

# ===================== conda environment =====================

source ~/.bashrc
micromamba activate fb-196-ic-0318-oe-2022-crit

micromamba env export > fitting-environment.yaml

export OE_LICENSE="/home/lilywang/oe_license.txt"

mkdir working-directory

ForceBalance.py optimize.in > force_balance.log
