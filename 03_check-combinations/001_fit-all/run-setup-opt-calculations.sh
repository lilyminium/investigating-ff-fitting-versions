#!/bin/bash
#
# Set the job name and wall time limit
#BSUB -J setup-opt-calculations
#BSUB -W 168:00
#
# Set the output and error output paths.
#BSUB -o  setup-opt-calculations-%J.o
#BSUB -e  setup-opt-calculations-%J.e
#
# Set any cpu options.
#BSUB -M 2

# ===================== conda environment =====================

source ~/.bashrc
micromamba activate fit-virtual-sites-tk010-py39

# python setup-opt-calculations.py                        \
#     --target-directory      targets                     \
#     --number                1                           \
#     --n-reps                3                           \
#     --environment           fb-193-tk-010-oe-2022       \
#     --environment           fb-193-tk-010-oe-2022-reordered     \
#     --environment           fb-195-tk-013-oe-2022-interchange-replace-cache     \
#     --environment           fb-196-ic-0318-oe-2022


# python setup-opt-calculations.py                        \
#     --target-directory      targets                     \
#     --number                1                           \
#     --n-reps                3                           \
#     --environment           fb-195-tk-013-oe-2022-interchange-replace-cache-switching     \

python setup-opt-calculations.py                        \
    --target-directory      targets                     \
    --number                1                           \
    --n-reps                3                           \
    --environment           fb-196-ic-0318-oe-2022-crit
