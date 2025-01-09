#! /bin/sh
#SBATCH --mem=5000M
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --job-name=genpot
date
/WORK/xyzhong_5/MATSv2/genpot3/genpot < genpot.in 
wait
date
