#!/bin/bash
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --nodes=1
#SBATCH --job-name=mult5-1
#SBATCH --exclusive 

echo "Starting" `date`

# script starts here


qxmin=-20.0
qxmax=20.0
qxstp=5.0

qymin=-20.0
qymax=20.0
qystp=5.0

icount=0
for i in $(seq $qxmin $qxstp $qxmax)
do
  for j in $(seq $qymin $qystp $qymax)
  do

#  echo $i $j 
    icount=$((icount+1))
    rm -rf $icount
  done
done

echo "Finished" `date`
