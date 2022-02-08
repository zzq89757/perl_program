use strict;
my $fa = shift;
my %hash;
my ($count_poly,$count_repeat,$count_mutation,$count_all);
open IN,"$fa";
while(my $id1 = <IN>){
  my $seq1 = <IN>;
  my $id2 = <IN>;
  my $seq2 = <IN>;
  my $rid = (split/_/,$id1)[0];
  $hash{data}{$rid}="$id1$seq1$id2$seq2";
  $count_all++;
  #寻找poly结构 此处划分标准为四个及以上相同碱基
  if($seq1=~/([AGCT])\1{3}/ ||$seq2=~/([AGCT])\1{3}/){
    # print "$id1$seq1";
    $count_poly++;
    next;
  }
  # if($seq1=~/([ACGT][ACTG])\1{2}/ || $seq2=~/([ACGT][ACTG])\1{2}/){
  if($seq1=~/([ACGT]{2,30})\1{2}/ || $seq2=~/([ACGT]{2,30})\1{2}/){
    $count_repeat++;
        # print "$id1$seq1";
    next;
  }
  my $msg1 = (split/_/,$id1)[-1];
  my $msg2 = (split/_/,$id2)[-1];
  my @base1 = $msg1 =~ /\d+([AGCT])/g;
  for my $item (@base1){
    if(!exists $hash{base1}{$item}){
      $hash{base1}{$item}=1
    }
  }
  my @sort1 = sort keys %{$hash{base1}};
  delete $hash{base1};
  my @base2 = $msg2 =~ /\d+([AGCT])/g;
  for my $item (@base2){
    if(!exists $hash{base2}{$item}){
      $hash{base2}{$item}=1
    }
  }
  my @sort2 =  sort keys %{$hash{base2}};
  delete $hash{base2};
  #此处判断存在漏检
  if(@base1>1&&@sort1==1 || @base2>1&&@sort2==1){
    # print "$id1$seq1";
    $count_mutation++;
  }
  else{
    print STDERR $hash{data}{$rid}
  }
}
my $pct_poly = sprintf "%.2f%%",100 * $count_poly/$count_all;
my $pct_re = sprintf "%.2f%%",100 * $count_repeat/$count_all;
my $pct_mu = sprintf "%.2f%%",100 * $count_mutation/$count_all;
my $ot = $count_all-$count_poly-$count_repeat-$count_mutation;
my $pct_ot = sprintf "%.2f%%",100* $ot/$count_all;
print "Type\tAll_Insert\tPoly\tTandem_Repeat\tRegular_Mutation\tUnknow\nCount\t$count_all\t$count_poly($pct_poly)\t$count_repeat($pct_re)\t$count_mutation($pct_mu)\t$ot($pct_ot)\n"