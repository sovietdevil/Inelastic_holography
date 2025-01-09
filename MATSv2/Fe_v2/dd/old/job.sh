#!/bin/bash

echo "Starting" `date`
# script starts here
#mkdir thick_files


  line=$( cat shifts | head -n 1 | tail -n 1 )
  read ox oy <<< $line
  echo $ox $oy

  echo $ox $oy > fort.75
  echo 48 48 >> fort.75

  for i in thick_files/* ; do
    echo $i
    j=$( basename $i )
    echo $j
  
    cp -f $i Fe.multislice
  
    cat insel2_template_thick | sed -e s/TTTT/${j}/ > insel2_template
#    if [ ! -f ../dd_res/${l}_e2_${j} ]; then
#      ./job_master3.pl
#      mv e1 ../dd_res/${l}_e1_${j}
#      mv e2 ../dd_res/${l}_e2_${j}
#    fi

  done # for thick_files


echo "Finished" `date`
