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
my $countFile = shift;
my $head = `samtools view -H $bamFile`;
my $totalReads = 0;
my $count_multi = 0;
my $count_unmap = 0;
my $count_abnormal = 0;
my $count_normal = 0;
my $max = 808;
# my $max = 933;
my $hash;
open COUNT,">$countFile";
# open INSERT,">./$bamFile.INSERT.sam";
# print INSERT $head;
# open CONTIG,">./$bamFile.CONTIG.sam";
# print CONTIG $head;
# open SA,">./$bamFile.SA.sam";
# print SA $head;
# open FR,">./$bamFile.FR.sam";
# print FR $head;
# open OTHERS,">./$bamFile.OTHERS.sam" or die "fuck\n";
# print OTHERS $head;
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to sure that samtools is available \n";
my $r1Chimeras;
my $r2Chimeras;
my $chimerasDenominator;
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
my $count_insert;
my $count_contig;
my $count_satag;
my $count_notfr;
my $count_others = 0;
my $Total;
my %hash;
while(<BAM>){
    chomp;
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
    my $getSecondaryFlag = $binFlag[8];
    my $getReadFailsVendorQualityCheckFlag = $binFlag[9];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    my $getReadMappedFlag = $getReadUnmappedFlag == 0? 1:0;
    my $getMateMappedFlag = $getMateUnmappedFlag == 0? 1:0;
    # next if(!$getReadPairedFlag || $binFlag[8] || $getSupplementaryAlignmentFlag || $getReadFailsVendorQualityCheckFlag || $getReadUnmappedFlag);
    my $negativeStart;
    my $positiveStart;
    my $matchLength;
    my $distance;
    $hash{$line[0]}{data}.="$_\n";

        #判断方向
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
    ++ $Total;
    if(!exists $hash{$line[0]}{insert}){
        if(abs $insertSize > $max){
            $hash{$line[0]}{insert} = 1;
        }
    }else{
        ++ $count_insert;
        # print INSERT $hash{$line[0]}{data};
    }
    if(!exists $hash{$line[0]}{contig}){
        if($line[6]!~/=/){
            $hash{$line[0]}{contig} = 1;
        }
    }else{
        ++ $count_contig;
        # print CONTIG $hash{$line[0]}{data};
    }
    if(!exists $hash{$line[0]}{satag}){
        if($line[11]=~/SA:Z:/){
            $hash{$line[0]}{satag} = 1;
        }else{
            $hash{$line[0]}{satag} = 0
        }
    }elsif ($hash{$line[0]}{satag}==1){
        ++ $count_satag;
        # print SA $hash{$line[0]}{data};
    }else{
        if($line[11]=~/SA:Z:/){
            ++ $count_satag;
            $hash{$line[0]}{satag} = 1;
            # print SA $hash{$line[0]}{data};
        }
    }
    if(!exists $hash{$line[0]}{FR}){
        if($distance==0){
            $hash{$line[0]}{FR} = 1;
        }else{
            $hash{$line[0]}{FR} = 0
        }
    }elsif ($hash{$line[0]}{FR}==1){
            ++ $count_notfr;
            # print FR $hash{$line[0]}{data};
    }else{
        if($distance==0){
        ++ $count_notfr;
        $hash{$line[0]}{FR} = 1;
        # print FR $hash{$line[0]}{data};
        }
    }
    if(!exists $hash{$line[0]}{not4}){
        if( $hash{$line[0]}{FR} ==0 && $hash{$line[0]}{satag} == 0 && !$hash{$line[0]}{contig} && !$hash{$line[0]}{insert}){
            $hash{$line[0]}{not4} = 1;
        }else{
            $hash{$line[0]}{not4} = 0
        }
    }elsif ($hash{$line[0]}{not4} == 1){
        if($hash{$line[0]}{FR} ==0 && $hash{$line[0]}{satag} == 0 && !$hash{$line[0]}{contig} && !$hash{$line[0]}{insert}){
        ++ $count_others;
        # print OTHERS $hash{$line[0]}{data};
        }
    }

}
my $total_insert = $Total/2;
print COUNT "TOTAL_ABNORMAL_INSERT\tUnexpctedInsertSize\tDifferentContig\tWithSATag\tUnexpectedOrientation\tOTHERS\n";
print COUNT "$total_insert\t$count_insert\t$count_contig\t$count_satag\t$count_notfr\t$count_others"