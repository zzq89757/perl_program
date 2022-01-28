use strict;
my $ref = shift;
my $chr = shift;
my $pos = shift;
my $cut_length = shift;
open IN,"$ref";
$/ = ">";
while(<IN>){
  chomp;
  my $ID = (split(/\n/,$_,2))[0];
  my @line = split/\s+/,$ID;
  next if($ID !~/\S+/ ||$line[0]!=$chr);
  if($line[0] == $chr){
    my $seq = (split(/\n/,$_,2))[1];
    $seq=~s/\n//g;
    my $cut = substr($seq,$pos-1,$cut_length);
    print $cut;
    last;
  }
}