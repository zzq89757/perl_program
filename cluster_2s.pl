use strict;
my $fa = shift;
# my 
open IN,"$fa";
while(my $msg = <IN>){
  chomp;
  my $seq = <IN>;
  my @line = split/_/,$msg;
  my $name = $line[0];
  my $cigar = $line[1];

  print $line[0];
}