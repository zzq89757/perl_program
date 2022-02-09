use strict;
my $blast = shift;
my $fasta = shift;
open IN,"$blast";
while(<IN>){
  chomp;
  my @line = split/\t/;
  my $msg = $line[0];
  my $score = $line[2];
  my $evalue = $line[10];
  my $chr = (split/_/,$msg)[2];
  my $pos = (split/_/,$msg)[3];
  my $sub_start = $line[8] < $line[9] ? $line[8]:$line[9];
  my $sub_end = $line[9] > $line[8] ? $line[9]:$line[8];
  my $oriente = $line[9] > $line[8] ? "+":"-";
  my $sub_length = $sub_end - $sub_start + 1 ;
  my $raw_cut = `perl /mnt/data/Users/zzq/Abnormal_align/zzq/script/sub_ref.pl /mnt/data/Users/zzq/Abnormal_align/zzq/somatic-b37_Homo_sapiens_assembly19.fasta $chr $sub_start $sub_length`;
  my $hubu_cut = $raw_cut =~ tr/atcgATCG/tagcTAGC/;
  my $cut = $oriente=="+" ?  $raw_cut:scalar reverse $hubu_cut;
  my $raw_fa = `perl /mnt/data/Users/zzq/Abnormal_align/zzq/script/pull_by_id.pl $fasta $msg`;
  # print "$msg\n$raw\n$sub_start\t$sub_end\n$cut\n";
  print ">$msg\n$raw_fa\n$line[8]\t$line[9]\t$oriente\t$score\t$evalue\n$cut\n";
}