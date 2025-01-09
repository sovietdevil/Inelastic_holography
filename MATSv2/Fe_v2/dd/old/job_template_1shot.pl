#!/usr/bin/perl

#chdir('WORKDIR');

$HOME="/home/rusz";

$taskfile = $ARGV[0];

$jobno = JOBNO;
$case = "CASE";
$scratch = "SCRATCH";

$ddbin="DDBIN";
$mdffbin="MDFFBIN";

open state, ">state_$jobno";

system("cp","${case}.multislice","${case}.multislice_${jobno}");
system("cp","insel2_template","insel2_template_${jobno}");
system("cp","${case}.struct","${case}.struct_${jobno}");
system("mkdir","${jobno}");

print state `date`;
open task, "<$taskfile";
while($taskline=<task>) {
  ($ii, $jj) = split /\s+/, $taskline;
  print state "$ii $jj     ...     dyndif started:"; print state `date`;
  open temp, "<insel2_template_${jobno}";
  open script, ">${case}.insel2_${jobno}";
  while($line=<temp>) {
    $line =~ s/XXXX/$ii/;
    $line =~ s/YYYY/$jj/;
    print script $line;
  }
  close script;
  close temp;
  `./$ddbin dyndif.def_$jobno > ddout_${ii}_${jj}`;
  rename "${jobno}/${case}.ddscstot2_$jobno", "${jobno}/${case}.ddscstot2_${ii}_${jj}";
}
print state "ALL DONE! "; print state `date`;
close state;

#`tar -cvzf ddscs_${jobno}.tar.gz ${jobno}/* ; cp ddscs_${jobno}.tar.gz . ; tar -xvf ddscs_${jobno}.tar.gz `;
`cp ${jobno}/* .`;
