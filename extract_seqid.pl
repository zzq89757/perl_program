use strict;
my $input = shift;
my $out = shift;
open OU,">>$out";
open IN,"$input";
while(<IN>){
  chomp;
  if(/^Seq/){
    my $seqid = (split/\s+/)[1];
    print OU "$seqid\n";
  }
}