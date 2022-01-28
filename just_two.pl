use strict;
use Statistics::Descriptive;
use List::Util qw/sum/;
use Getopt::Long;
my $bamFile = shift;
my $max = shift;
my $abnormal = shift;
my $r1Reads;
my $r2Reads;
my $Total;
my $PairReads;
my $diffContig;
my $abnormalInsert;
my $orient;
my %hash;
my @insertSize;
my $normal;
my $countnormal;
my $count_abnormal;
my $cleanzero;
my $count_zero;
open IN,"samtools view $bamFile|";
open ABNORM,">$abnormal";
while(<IN>){
    chomp;
    my @line = split /\s+/;
    my $flag = $line[1];
    my $position = $line[3];
    my $cigar = $line[5];
    my $matePosition = $line[7];
    my $insertSize = $line[8];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag = split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getReadPairedFlag = $binFlag[0];
    my $getReadUnmappedFlag = $binFlag[2];
    my $getMateUnmappedFlag = $binFlag[3];
    my $getReadReverseFlag = $binFlag[4];
    my $getMateReverseFlag = $binFlag[5];
    my $getFirstPairFlag = $binFlag[6];
    my $getSecondaryFlag = $binFlag[8];
    my $getReadFailsVendorQualityCheckFlag = $binFlag[9];
    my $getDupFlag = $binFlag[10];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    my $getReadMappedFlag = $getReadUnmappedFlag == 0? 1:0;
    my $getMateMappedFlag = $getMateUnmappedFlag == 0? 1:0;
    if($getFirstPairFlag){
        $r1Reads++;
    }else{
        $r2Reads++
    }
        my $positiveStart;
        my $negativeStart;
        my $orientation;
        #get tandem
        if($getReadReverseFlag == $getMateReverseFlag){
           $orientation = "TANDEM";
        }else{
            if($getReadReverseFlag){
                my $matchLength = sum ($cigar =~/(\d+)[DMNX]/g);
                $positiveStart = $matePosition;
                $negativeStart = $position + $matchLength -1;
            }else{
                $positiveStart = $position;
                $negativeStart = $position + $insertSize;
            }
            $orientation = $negativeStart > $positiveStart ? "FR" : "RF";
        }
    # $rawzero .= "$_\n" if($insertSize == 0);
    $normal .= "$_\n" if($getReadMappedFlag && $getMateMappedFlag && $cigar!~/S/ && abs $insertSize <= $max);
    $countnormal++ if($getReadMappedFlag && $getMateMappedFlag && $cigar!~/S/ && abs $insertSize <= $max);
    next if($getReadMappedFlag && $getMateMappedFlag && $cigar!~/S/ && abs $insertSize <= $max);
    $count_abnormal++;
    $cleanzero .= "$_\n" if($insertSize == 0);
    $count_zero++ if($insertSize == 0);
    # $abnormal.= "$_\n";
    print ABNORM "$_\n";
}