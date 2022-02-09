use strict;
my $fa = shift;
my %hash;
open IN,"$fa";
while(my $id = <IN>){
  my $seq = <IN>;
  my @msg = split/_/,$id;
  my $base = $msg[-1];
  next if($base!~/[AGCT]/ ||$seq=~/([AGCT])\1{3}/);
  my @bases = $base =~ /\d+([AGCT])/g;
  for my $item (@bases){
    if(!exists $hash{base}{$item}){
    $hash{base}{$item}=1
    }
  }
  my @sort =  sort keys %{$hash{base}};
  delete $hash{base};
  next if (@sort > 2||@sort ==0);
  if(@bases>1&&@sort==1){
    print "$id";
  }
}