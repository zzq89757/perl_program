use strict;
use List::Util qw(min max);
my ($count_nos,$total,$count_s,$pct_s,$pct_nos,%hash);
open BAM,"samtools view P6_1_12178_ND608_N616_2.sort.markdup.bam.OTHERS.sam|";
while(<BAM>){
    chomp;
    my @line = split /\s+/;
    my $cigar = $line[5];
    $total++;

    if($cigar=~/(\d+)S/){
      my @s_num = $cigar=~/(\d+)S/g;
      my $max = max @s_num;
      $count_s++;
      $hash{length}{$max}++;
    }else{
      $count_nos++;
    }
}
$pct_s=sprintf "%.2f%%",100 * $count_s/$total;
$pct_nos=sprintf "%.2f%%",100 * $count_nos/$total;
print "sample\ttotal\twithS\twithoutS\n";
print "P6-1\t$total\t$count_s($pct_s)\t$count_nos($pct_nos)\n";
foreach my $key (sort{$hash{length}{$b} <=> $hash{length}{$a}} keys %{$hash{length}}){
  print "$key\t$hash{length}{$key}\n";
}