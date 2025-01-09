#!/usr/bin/perl

use POSIX;


# shall we split the real/imag parts?

$qxmin = -20.0;
$qxmax =  20.0;
$qxstp =  5.0;

$qymin = -20.0;
$qymax =  20.0;
$qystp =  5.0;

# extract results

open res, ">e1";
for($i=$qxmin;$i<=$qxmax;$i+=$qxstp) {
  for($j=$qymin;$j<=$qymax;$j+=$qystp) {
    $fname = sprintf "Fe.ddscstot2_%-3.1f_%-3.1f", $i, $j, $op;
#    print "$fname\n";
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

