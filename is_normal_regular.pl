use strict;
my $bam = shift;
open IN,"samtools view $bam|";
while(<IN>){
  /MD:Z:(\S+)/;
  # print $1;
  my @c = $1=~/\d+([AGCT])/g;
  # print @c;
  next if(@c<3);
  print;
  # my @line = split/\t/;
}