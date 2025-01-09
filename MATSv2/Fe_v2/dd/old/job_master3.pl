#!/usr/bin/perl

use POSIX;

$workdir = `pwd`;
chomp $workdir;
print "$workdir\n";
#chdir($workdir);

# shall we split the real/imag parts?
$reim=0;

#print "|$jobcount|\n";
$case = "Fe";

#$ddbin="dyndif_v15_fe_sf";
$ddbin="dyndif_v15_fe_sf_multbw_opmaps6_zloc";
$mdffbin="mdff_dip";

$qxmin = -20.0;
$qxmax =  20.01;
$qxstp =  4.0;

$qymin = -20.0;
$qymax =  20.01;
$qystp =  4.0;

# read ".machines" file
open mach, "<.machines";
$l = <mach>;
$jobcount = 0;
while($l =~ m/^1:(.*)$/) {
  $jobcount++;
  $machine[$jobcount] = $1;
  print "$machine[$jobcount]\n";
  $l = <mach>;
}
close mach;

print $jobcount;

# create definition files
for($i=0;$i<$jobcount;$i++) {
  open def, ">dyndif.def_$i";
  print def "16,'${case}.insel2_${i}',   'old',    'formatted',0\n";
#  print def "17,'${case}.glist',    'unknown','formatted',0\n";
#  print def "18,'${case}.glist2',   'unknown','formatted',0\n";
  print def "19,'${case}.glist3_${i}',   'unknown','formatted',0\n";
  print def "20,'${case}.struct_${i}',       'old',    'formatted',0\n";
#  print def "25,'${case}.w2kpot',        'unknown' 'formatted',0\n";
#  print def "40,'${case}.qqprlist_${i}', 'unknown','formatted',0\n";
#  print def "60,'${case}.mdff_${i}',     'unknown','formatted',0\n";
  print def "74,'${case}.multislice_${i}','unknown','formatted',0\n";
#  print def "80,'${case}.ddscstot_${i}', 'unknown','formatted',0\n";
  print def "82,'${i}/${case}.ddscstot2_${i}','unknown','formatted',0\n";
#  print def "50,'${case}.blochs_${i}',   'unknown','binary',   0\n";
  close def;
}
  
#for($i=0;$i<$jobcount;$i++) {
#  open def, ">mdff.def_$i";
#  print def "15,'${case}.qqprlist_${i}', 'old',    'formatted',0\n";
#  print def "16,'${case}.insel2_${i}',   'old',    'formatted',0\n";
#  print def "37,'${case}.inmdffdip_${i}','old',    'formatted',0\n";
#  print def "60,'${case}.mdff_${i}',     'unknown','formatted',0\n";
#  close def;
#}

# create tasklist files
for($i=0;$i<$jobcount;$i++) {
  open $tmp, '>', "tasklist_${i}";
  push @task, $tmp;
  print "$mp\n";
  undef $tmp;
}
$jobno=0;
for($i=$qxmin;$i<$qxmax;$i+=$qxstp) {
  for($j=$qymin;$j<$qymax;$j+=$qystp) {
    $jobno = $jobno % $jobcount;
#    print "$jobno\n";
    printf {$task[$jobno]} "%05.2f %05.2f\n", $i, $j;
    $jobno++;
  }
}
for($i=0;$i<$jobcount;$i++) {
  close $task[$i];
}

# create individual job scripts
if($sx) {
  $stp = $qxstp;
  $min = $qxmin;
} else {
  $stp = $qystp;
  $min = $qymin;
}

for($i=0;$i<$jobcount;$i++) {
  open temp, "<job_template_1shot.pl";
  open script, ">job_${i}.pl";
  while($line=<temp>) {
    $line =~ s/WORKDIR/$workdir/;
    $line =~ s/JOBNO/$i/;
    $line =~ s/CASE/$case/;
    $line =~ s/DDBIN/$ddbin/;
    $line =~ s/MDFFBIN/$mdffbin/;
    $line =~ s/SCRATCH/${ENV{'SNIC_TMP'}}/;
    print script $line;
  }
  close script;
  close temp;
  chmod 0744, "job_${i}.pl";
}

# run job scripts in parallel
open machines, "<.machines";
my @children;
for($i=0;$i<$jobcount;$i++) {
  $line = <machines>;
  chomp $line;
  @arr = split /:/, $line;
  print "$arr[1]\n";
  $command = "(cd $workdir; touch .lockfile_${i}; ./job_${i}.pl tasklist_${i}; rm -f .lockfile_${i})";
  print "$command\n";
  my $child = fork;
  exec "$command" if $child==0;
  push @children, $child;
  sleep 1;
}
close machines;

# wait until all jobs are finished
foreach $child (@children) { 
  waitpid $child, 0; 
}
#wait;

# extract results

open res, ">e1";
for($i=$qxmin;$i<$qxmax;$i+=$qxstp) {
  for($j=$qymin;$j<$qymax;$j+=$qystp) {
    $fname = sprintf "${case}.ddscstot2_%05.2f_%05.2f", $i, $j, $op;
    open inp, $fname;
    $s = "";
    $l = <inp>;
    while($l!~m/^\s*$/) { 
      $l = <inp>;
      @a = split /\s+/, $l;
      $s .= "   ".$a[2]."  ".$a[3]."  ".$a[4]."  ".$a[5]."  ".$a[6]."  ".$a[7]."  ".$a[8]."  ".$a[9]."  ".$a[10];
    }
    print res   "$s\n";
    close inp;
  }
}
close res;

open res, ">e2";
for($i=$qxmin;$i<$qxmax;$i+=$qxstp) {
  for($j=$qymin;$j<$qymax;$j+=$qystp) {
    $fname = sprintf "${case}.ddscstot2_%05.2f_%05.2f", $i, $j, $op;
    open inp, $fname;
    $s = "";
    $l = <inp>;
    while($l!~m/^\s*$/) { 
      $l = <inp>; 
    }
    $l = <inp>;
    while($l!~m/^\s*$/) { 
      $l = <inp>;
      @a = split /\s+/, $l;
      $s .= "   ".$a[2]."  ".$a[3]."  ".$a[4]."  ".$a[5]."  ".$a[6]."  ".$a[7]."  ".$a[8]."  ".$a[9]."  ".$a[10];
    }
    print res   "$s\n";
    close inp;
  }
}
close res;
