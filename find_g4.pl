use strict;
my $bam = shift;
my $cut = shift;
my %hash;
open IN,"samtools view $bam|";
while(<IN>){
  chomp;
  my @line = split;
  my $name = $line[0];
  my $chr = $line[2];
  my $position = $line[3];
  my $cigar = $line[5];
  my $seq = $line[9];
  my $mdtag = "----";
  $mdtag = $1 if(/MD:Z:(\S+)/);
  my @array = $cigar =~/(\d+)S/g;
  if(!exists $hash{S}{$name}){
    if((@array == 1 && $array[0]>4) || (@array ==2 && ($array[0]>4 || $array[1]>4)) ){
      $hash{S}{$name} = 1;
      $hash{data}{$name} .= ">$name"."_$cigar"."_$chr"."_$position"."_$mdtag"."\n$seq\n";
    }else{
      $hash{data}{$name} .= ">$name"."_$cigar"."_$chr"."_$position"."_$mdtag"."\n$seq\n";
      $hash{S}{$name} = 0;
    }
  }elsif($hash{S}{$name}==1){
    $hash{data}{$name} .= ">$name"."_$cigar"."_$chr"."_$position"."_$mdtag"."\n$seq\n";
    delete $hash{S}{$name};
  }elsif($hash{S}{$name}==0){
    if((@array == 1 && $array[0]>4) || (@array ==2 && ($array[0]>4 || $array[1]>4)) ){
      $hash{data}{$name} .= ">$name"."_$cigar"."_$chr"."_$position"."_$mdtag"."\n$seq\n";
    }else{
      delete $hash{data}{$name}
    }
    delete $hash{S}{$name}
  }
}
delete $hash{S};
for my $key(keys %{$hash{data}}){
  print $hash{data}{$key}
}