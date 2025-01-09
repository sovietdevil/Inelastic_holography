#!/usr/bin/perl

sub removefiles {
  my $nj = shift;
  my $fn = shift;
  my $sf = shift;
  $str = '';
  for($i=0;$i<$nj;$i++) {
    $str .= " ${fn}_${i}${sf}";
  }
  print "Removing files ${fn}_*${sf} ...\n";
#  print "$str\n";
  print `rm -f $str`;
}

# read ".machines" file
open mach, "<.machines";
$l = <mach>;
$jobcount = 0;
while($l =~ m/^1:(.*)$/) {
  $jobcount++;
  $machine[$jobcount] = $1;
#  print "$machine[$jobcount]\n";
  $l = <mach>;
}
close mach;

print "Number of parallel CPUs: $jobcount\n";

#$case="hcpCo_orth";
#$case="Fe";
#$case="afm_LaOMnAs";
#$case="LaMnO3";
#$case="FePt";
$case="Fe";

`rm -f ${case}.ddscstot*`;
`rm -f ddout_*`;
#`rm fort.*`;
`rm ${case}.glist*`;
`rm ${case}.w2kpot`;
`rm insel2_template`;

# remove files
&removefiles($jobcount,"dyndif.def","");
&removefiles($jobcount,"mdff.def","");
&removefiles($jobcount,"tasklist","");
&removefiles($jobcount,"state","");
&removefiles($jobcount,"*.mdff","");
&removefiles($jobcount,"*.mdffim","");
&removefiles($jobcount,"*.mdffre","");
&removefiles($jobcount,"*.inmdffdip","");
&removefiles($jobcount,"*.outputmdff","");
&removefiles($jobcount,"*.qqprlist","");
&removefiles($jobcount,"*.insel2","");
&removefiles($jobcount,"*.glist","");
&removefiles($jobcount,"*.glist2","");
&removefiles($jobcount,"*.blochs","");
&removefiles($jobcount,"*.multislice","");
&removefiles($jobcount,"job",".pl");
&removefiles($jobcount,"ddout","");
&removefiles($jobcount,".lockfile","");
