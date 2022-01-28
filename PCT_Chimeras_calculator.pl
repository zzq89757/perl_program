# -------------------- about binFlag -----------------------------
# 0:read paired
# 1:read mapped in proper pair  ****
# 2:read unmapped
# 3:mate unmapped
# 4:read reverse strand
# 5:mate reverse strand
# 6:first in pair
# 7:second in pair
# 8：not primary alignment
# 9:read fails platform/vendor quality checks
# 10:read is PCR or optical duplicate  ****
# 11:supplementary alignment
use strict;
use Getopt::Long;
sub usage {
        print STDERR << "USAGE";
Description     PCT_Chimeras_calculator By zzq
Author                  zhengzhiqiang\@vazyme.com
Data                    2021.11
Version                 V1.0
Usage                   perl $0 bamFile


 e.g.:
         perl $0  /xxx/xxx.bam
USAGE

exit 0
}
usage() if(@ARGV<1);
my $bamFile = shift;
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to make sure that samtools is available \n";
my $r1Chimeras;
my $r2Chimeras;
my $chimerasDenominator;
my $r1Reads;
my $r2Reads;
my $totalReads;
my $r1ChimerasDenominator;
my $r2ChimerasDenominator;
my $r1UnexpctedDistance;
my $r1DifferentContig;
my $r1WithSATag;
my $r1UnexpectedOrientation;
my $r2UnexpctedDistance;
my $r2DifferentContig;
my $r2WithSATag;
my $r2UnexpectedOrientation;
my %hash;
while(<BAM>){
    my @line = split /\t/;
    my $flag = $line[1];
    my $position = $line[3];
    my $mappingQuality = $line[4];
    my $cigar = $line[5];
    my $matePosition = $line[7];
    my $insertSize = $line[8];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getReadPairedFlag = $binFlag[0];
    my $getReadUnmappedFlag = $binFlag[2];
    my $getMateUnmappedFlag = $binFlag[3];
    my $getReadReverseFlag = $binFlag[4];
    my $getFirstPairFlag = $binFlag[6];
    my $getReadFailsVendorQualityCheckFlag = $binFlag[9];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    #bin[8]和$getSupplementaryAlignmentFlag不通过则总reads数不增加 最后两个指标无论是否通过都加总reads数
        ++$totalReads;
    next if($binFlag[8] || $getSupplementaryAlignmentFlag);
    next if($getReadFailsVendorQualityCheckFlag || $getReadUnmappedFlag || !$getReadPairedFlag);
	if(!$getMateUnmappedFlag){
        ++ $chimerasDenominator;
        my $positiveStart;
        my $negativeStart;
        my $matchLength;
        my $distance;
        if($binFlag[4] == $binFlag[5]){
            $distance = 0;
        }else{
            if($getReadReverseFlag){
                my @match = $line[5] =~/(\d+)[DMNX]/g;
                print STDERR "$line[1] cigar erro!!\n" if(@match eq 0);
                foreach my $i(@match){
                    $matchLength += $i;
                }
                $positiveStart = $matePosition;
                $negativeStart = $position + $matchLength -1;
            }else{
                $positiveStart = $position;
                $negativeStart = $position + $insertSize;
            }
            $distance = $negativeStart > $positiveStart ? 1 : 0;
        }
        ++ $r1ChimerasDenominator if($getFirstPairFlag);
        ++ $r2ChimerasDenominator if(!$getFirstPairFlag);        
            if($line[6]!~/=/ || abs $insertSize >100000 || $line[11]=~/SA:Z:/ || !$distance){
                if($getFirstPairFlag){
                    ++ $r1Chimeras;
                    ++ $r1UnexpctedDistance if(abs $insertSize >100000);
                    ++ $r1UnexpectedOrientation if(!$distance);
                    ++ $r1DifferentContig if($line[6]!~/=/);
                    ++ $r1WithSATag if($line[11] =~/SA:Z:/);
                }else{           
                    ++ $r2Chimeras;
                    ++ $r2UnexpctedDistance if(abs $insertSize >100000);
                    ++ $r2UnexpectedOrientation if(!$distance);
                    ++ $r2DifferentContig if($line[6]!~/=/);
                    ++ $r2WithSATag if($line[11] =~/SA:Z:/);
                }
            }
    }else{
        if ($mappingQuality >= 20) {
            ++ $chimerasDenominator;
            ++ $r1ChimerasDenominator if($getFirstPairFlag);
            ++ $r2ChimerasDenominator if(!$getFirstPairFlag);
            if($line[11]=~/SA:Z:/) {
                if($getFirstPairFlag){
                    ++ $r1Chimeras;
                    ++ $r1WithSATag;
                }else{
                    ++ $r2Chimeras;
                    ++ $r2WithSATag;   
                }
            }
        }
    }
}
my $totalChimeras = $r1Chimeras + $r2Chimeras;
my $pctR1Chimeras  = sprintf "%.6f",$r1Chimeras/$r1ChimerasDenominator;
my $pctR2Chimeras  = sprintf "%.6f",$r2Chimeras/$r2ChimerasDenominator;
my $pctChimeras = sprintf "%.6f",($r1Chimeras+$r2Chimeras)/$chimerasDenominator;
print "-----------------Total Count-------------------\n";
print "COUNT_R1\t$r1Reads\nCOUNT_R2\t$r2Reads\nTOTAL_READS\t$totalReads\n";
print "-----------------Chimeras Count-------------------\n";
print "COUNT_CHIMERAS_OF_R1\t$r1Chimeras\nCOUNT_CHIMERAS_OF_R2\t$r2Chimeras\nCOUNT_CHIMERAS_OF_ALL\t$totalChimeras\nPCT_CHIMERAS_OF_R1\t$pctR1Chimeras\nPCT_CHIMERAS_OF_R2\t$pctR2Chimeras\nPCT_CHIMERAS_OF_ALL\t$pctChimeras\n";
print "---------------All Chimeras Type Count------------------\n";
print "ReadCount\\Type\tUnexpctedDistance\tDifferentContig\tWithSATag\tUnexpectedOrientation\n";
print "FirstOfPair\t$r1UnexpctedDistance\t$r1DifferentContig\t$r1WithSATag\t$r1UnexpectedOrientation\n";
print "SecondOfPair\t$r2UnexpctedDistance\t$r2DifferentContig\t$r2WithSATag\t$r2UnexpectedOrientation\n";