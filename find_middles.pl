use strict;
use List::Util qw(min max);
use List::MoreUtils qw(indexes);
my $bam = shift;
my ($count_nos,$total,$count_s,$pct_s,$pct_nos,%hash);
open BAM,"samtools view $bam|";
while(<BAM>){
    chomp;
    my @line = split /\s+/;
    my $cigar = $line[5];
    $total++;

    if($cigar=~/(\d+)S/){
      my @array = $cigar=~/[MDNXS]/g;
      my @indexes = indexes { $_ eq 'S' } @array;
      my @s_num = $cigar=~/(\d+)S/g;
        my $endix = @array -1;
        if(@s_num >2){
          print "ass"
        }
      for my $ind(@indexes){
        if($ind!=0 && $ind != $endix){
        print "@indexes\t$ind\t$endix";
        }
      }
      # print "\n";
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
  # print "$key\t$hash{length}{$key}\n";
}