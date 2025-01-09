#!/bin/bash
#SBATCH --mem=64000M
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --job-name=mult5-1
#SBATCH --exclusive 

date
/WORK/xyzhong_5/MATSv2/mult6_tilt/mult6 < mult6.in 
date
