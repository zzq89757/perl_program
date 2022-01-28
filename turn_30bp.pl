use strict;
use List::Util qw(min max);
use List::MoreUtils qw(indexes);

my $bam = shift;
my @prefix = split/\_/,$bam;
my $prefix = $prefix[0]."_".$prefix[1];
my ($count_nos,$total,$count_s,$pct_s,$pct_nos,%hash);
open BAM,"samtools view $bam|";
open FA,">$prefix.fasta";
open COUNT,">$prefix.count";
while(<BAM>){
    chomp;
    my @line = split /\s+/;
    my $binnaryFlagString = sprintf("%b",$line[1]);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getFirstPairFlag = $binFlag[6];
    my $readname = $getFirstPairFlag ? $line[0]:$line[0]."_2";
    my $cigar = $line[5];
    my $seq = $line[9];
    my $qual = $line[10];
    $total++;
    if($cigar=~/(\d+)S/){
      my @array = $cigar=~/[MDNXS]/g;
      my @indexes = indexes { $_ eq 'S' } @array;
      my @s_num = $cigar=~/(\d+)S/g;
      my $max = max @s_num;
      my @max_index = indexes { $_ eq $max } @s_num;
      my $cut;
      #获取S最大长度的索引 若只有一个 判断array中S索引 若有2个 若max为0 则S在前 反之 则在最后
        my $slength = $s_num[0];
      if(@s_num ==1){
          if($indexes[0] ==0){
            #说明S在前面 
            #大于30则直接输出 否则进行整条提取
            $cut = $slength>30? $seq:substr($seq,0,30);
          }else{
            $cut = $slength>30? $seq:substr($seq,-30,30);
          }
      }else{
        if($max_index[0] ==0){
            $cut = $slength>30? $seq:substr($seq,0,30);
        }else{
          $slength = $s_num[1];
            $cut = $slength>30? $seq:substr($seq,-30,30);
        }
      }
      $count_s++;
      $hash{length}{$max}++;
      print FA ">$readname\n$cut\n";
    }else{
      $count_nos++;
    }
}
$pct_s=sprintf "%.2f%%",100 * $count_s/$total;
$pct_nos=sprintf "%.2f%%",100 * $count_nos/$total;
print COUNT "sample\ttotal\twithS\twithoutS\n";
print COUNT "P6-1\t$total\t$count_s($pct_s)\t$count_nos($pct_nos)\n";
foreach my $key (sort{$hash{length}{$b} <=> $hash{length}{$a}} keys %{$hash{length}}){
  print COUNT "$key\t$hash{length}{$key}\n";
}