use strict;
my $bamFile = shift;
my $list = shift;
open OU,">$list";
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to make sure that samtools is available \n";
while(<BAM>){
  chomp;
  my @line = split/\s+/;
  my $read_id = $line[0];
  print OU "$line[0]\n";
}