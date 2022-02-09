use strict;
my $fa = shift;
open IN,"$fa";
while(my $id = <IN>){
  my $seq = <IN>;
  if($seq=~/([AGCT])\1{3}/){
    print "$id$seq";
    next;
  }else{
  }
}