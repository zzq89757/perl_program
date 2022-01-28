use strict;
my $bam = shift;
my %hash;
open BAM,"samtools view -F 256 $bam|";
while(<BAM>){
  chomp;
  my @line = split /\s+/;
  if(!exists $hash{XA}{$line[0]}){
    if(/XA/){
      $hash{XA}{$line[0]}=1;
    }
  }else{
    if(/XA/){
      $hash{XA}{$line[0]}=2;
    }
  }
}
foreach my $key (keys %{$hash{XA}}){
  print "$key\n" if($hash{XA}{$key}==1);
}