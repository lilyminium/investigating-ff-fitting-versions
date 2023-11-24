#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J "worker_1_195[1-8]"
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  worker_1_195-%J.o
#BSUB -e  worker_1_195-%J.e
#
# Set any cpu options.
#BSUB -M 8
#BSUB -n 8

# ===================== conda environment =====================

source ~/.bashrc
micromamba activate fb-195-tk-013-oe-2022-interchange-replace-cache

export OE_LICENSE="/home/lilywang/oe_license.txt"

DESTINATION="lt17"

for i in $(seq 1 8); do
    work_queue_worker --cores 1 -s working-directory/ --disk-threshold=0.002 --disk=3000 --memory-threshold=1000 --memory 1000 -b 20 -t 3600 "${DESTINATION}.hpc.private:22001" &
done
wait
