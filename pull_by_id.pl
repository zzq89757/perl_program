use strict;
my $fa = shift;
my $msg = shift;
open IN,"$fa";
while(my $id = <IN>){
  my $seq = <IN>;
  chomp $seq;
  if($id=~/$msg/){
    print $seq;
    last;
  }
}