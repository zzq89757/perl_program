use strict;
my $fa = shift;
my($total_base,$count_gc,@gc,%hash);
open IN,"$fa";
while(my $id1 = <IN>){
  chomp;
  my $seq1 = <IN>;
  my $id2 = <IN>;
  my $seq2 = <IN>;
  $total_base+=length$seq1;
  $total_base+=length$seq2;
  my @gc1 = $seq1=~/[AT]/g;
  my @gc2 = $seq2=~/[AT]/g;
  my $pct1 = @gc1/length $seq1; 
  my $pct2 = @gc2/length $seq2;
  push (@gc,$pct1);
  push (@gc,$pct2); 
  # $count_gc+=@gc1;
  # $count_gc+=@gc2;
}
for my $item (@gc){
  if($item<0.3){
    $hash{3}++
  }elsif($item<0.4){
    $hash{4}++
  }elsif($item<0.5){
    $hash{5}++
  }elsif($item<0.6){
    $hash{6}++
  }elsif($item<0.7){
    $hash{7}++
  }else{
    $hash{8}++
  }
  # print $item."\n";
}
foreach my $key (sort {$a<=>$b} keys %hash){
  my $n = $key -1;
  print "0.$n-0.$key\t$hash{$key}\n"
}
# my $pct_gc = sprintf "%.2f%%",100* $count_gc/$total_base;
# print "$count_gc\t$total_base\t$pct_gc\n";