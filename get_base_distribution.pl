use strict;
use List::Util qw(min max);
use List::MoreUtils qw(indexes);

my $r1 = shift;
my $r2 = shift;
my @prefix = split/\_/,$r1;
my $prefix = $prefix[0];
my ($count_nos,$total,$count_s,$pct_s,$pct_nos,%hash);
open R1,"$r1";
open R2,"$r2";
my $count_a = 0;
my $count_g = 0;
my $count_t = 0;
my $count_c = 0;
while(my $msg = <R1>){
    chomp $msg;
    my @line = split /_/,$msg;
    my $cigar = $line[3];
    next if ($cigar!~/^2S/ && $cigar!~/2S$/);
    my $seq = <R1>;
    chomp $seq;
    my $cut;
    if($cigar=~/^2S/){
      $cut = substr($seq,0,2);
    }else{
      $cut = substr($seq,-2,2);
    }
    my $g = () = $cut=~ /G/g;
    my $c = () = $cut=~ /C/g;
    my $t = () = $cut=~ /T/g;
    my $a = () = $cut=~ /A/g;
    $count_a+=$a if($a);
    $count_g+=$g if($g);
    $count_c+=$c if($c);
    $count_t+=$t if($t);
}
close R1;
my $count_all = $count_a + $count_c + $count_g + $count_t;
my $pct_a = sprintf "%.2f%%",100 * $count_a/$count_all;
my $pct_g = sprintf "%.2f%%",100 * $count_g/$count_all;
my $pct_c = sprintf "%.2f%%",100 * $count_c/$count_all;
my $pct_t = sprintf "%.2f%%",100 * $count_t/$count_all;
print  "R1\tALL:$count_all\tA:$count_a($pct_a)\tG:$count_g($pct_g)\tC:$count_c($pct_c)\tG:$count_g($pct_g)\n";

$count_a = 0;
$count_g = 0;
$count_t = 0;
$count_c = 0;
while(my $msg = <R2>){
    chomp $msg;
    my @line = split /_/,$msg;
    my $cigar = $line[3];
    next if ($cigar!~/^2S/ && $cigar!~/2S$/);
    my $seq = <R2>;
    chomp $seq;
    my $cut;
    if($cigar=~/^2S/){
      $cut = substr($seq,0,2);
    }else{
      $cut = substr($seq,-2,2);
    }
    my $g = () = $cut=~ /G/g;
    my $c = () = $cut=~ /C/g;
    my $t = () = $cut=~ /T/g;
    my $a = () = $cut=~ /A/g;
    $count_a+=$a if($a);
    $count_g+=$g if($g);
    $count_c+=$c if($c);
    $count_t+=$t if($t);
}
$count_all = $count_a + $count_c + $count_g + $count_t;
$pct_a = sprintf "%.2f%%",100 * $count_a/$count_all;
$pct_g = sprintf "%.2f%%",100 * $count_g/$count_all;
$pct_c = sprintf "%.2f%%",100 * $count_c/$count_all;
$pct_t = sprintf "%.2f%%",100 * $count_t/$count_all;
print  "R2\tALL:$count_all\tA:$count_a($pct_a)\tG:$count_g($pct_g)\tC:$count_c($pct_c)\tG:$count_g($pct_g)\n";