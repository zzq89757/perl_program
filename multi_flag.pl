use strict;
my $bam = shift;
open BAM,"samtools view $bam|";
while(<BAM>){
    chomp;
    my @line = split /\s+/;
    my $flag = $line[1];
    my $position = $line[3];
    my $mappingQuality = $line[4];
    my $cigar = $line[5];
    my $matePosition = $line[7];
    my $insertSize = $line[8];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getSecondaryFlag = $binFlag[8];
    my $getSupplementaryAlignmentFlag = $binFlag[11];

    # not primary alignment
    # print "have not primary alignment flag but without XA :";
    # if($getSecondaryFlag){
      if($_=~/XS/){
              print "$_\n";
      }

    # }
    # print "have supplementary alignment flag but without XA :";
    # if($binFlag[11]){
    #   next if(/XA/);
    #   print "$_\n";
    # }
}
