use strict;
use List::Util qw(min max);
use List::MoreUtils qw(indexes);

my $r1 = shift;
my @prefix = split/\./,$r1;
my $prefix = $prefix[0];
my ($count_nos,$total,$count_s,$pct_s,$pct_nos,%hash);
open R1,"$r1";
open CUT1,">$prefix.cut.fasta";
while(my $msg = <R1>){
    chomp $msg;
    my @line = split /_/,$msg;
    my $cigar = $line[1];
    my $seq = <R1>;
    chomp $seq;
    my $cut;
    if($cigar=~/(\d+)S/){
      my @array = $cigar=~/[MDNXS]/g;
      my @indexes = indexes { $_ eq 'S' } @array;
      my @s_num = $cigar=~/(\d+)S/g;
      my $max = max @s_num;
      my @max_index = indexes { $_ eq $max } @s_num;
      # my $cut;
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
    }else{
      $cut = $seq;
    }
    print CUT1 "$msg\n$cut\n"
}