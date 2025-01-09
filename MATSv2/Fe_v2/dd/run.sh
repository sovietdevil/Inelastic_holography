#!/bin/bash
#SBATCH --ntasks=9
#SBATCH --ntasks-per-node=9
#SBATCH --nodes=1
#SBATCH --job-name=mult5-1
#SBATCH --exclusive 

echo "Starting" `date`

# script starts here

for i in thick_files/* ; do
  j=$( basename $i )
  cp -f $i Fe.multislice
  cat insel2_template_thick | sed -e s/TTTT/${j}/ > insel2_template
done # for thick_files

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
    mkdir $icount
    cat insel2_template | sed -e s/XXXX/$i/ | sed -e s/YYYY/$j/ > Fe.insel2 
    cp dyndif.def Fe.insel2 Fe.struct Fe.multislice fort.75 $icount
    cd $icount
      /WORK/xyzhong_5/MATSv2/SRC_mult_bw_opmaps6_zloc/dyndif_v15_fe_sf_multbw_opmaps6_zloc dyndif.def > ddout &
      if (( $icount % $SLURM_NPROCS == 0 ))
      then
        wait
      fi 
    cd ../
  done
done

echo "Finished" `date`
