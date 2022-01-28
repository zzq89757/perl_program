use strict;
my $fa = shift;
open IN,"$fa";
while(my $id1 = <IN>){
  chomp;
  my $seq1 = <IN>;
  my $id2 = <IN>;
  my $seq2 = <IN>;
  if($id1 eq $id2){
    print "$id1$seq1$id2$seq2";
  }
}