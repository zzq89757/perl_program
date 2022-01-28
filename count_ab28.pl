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
Description     Count_abnormal By zzq
Author                  zhengzhiqiang\@vazyme.com
Data                    2021.11
Version                 V1.0
Usage                   perl $0 bamFile countFile


 e.g.:
         perl $0  /xxx/xxx.bam xxx.count
USAGE

exit 0
}
usage() if(@ARGV<2);
my $bamFile = shift;
my $countFile = shift;
my $head = `samtools view -H $bamFile`;
my $totalReads = 0;
my $count_xa = 0;
my $count_multi = 0;
my $count_unmap = 0;
my $count_abnormal = 0;
my $count_normal = 0;
my $max = 808;
# my $max = 933;
my $hash;
open COUNT,">$countFile";
open XA,">./$bamFile.XA.sam";
print XA $head;
open MULTI,">./$bamFile.MULTI.sam";
print MULTI $head;
open UMAP,">./$bamFile.UMAP.sam";
print UMAP $head;
open NORMAL,">./$bamFile.NORMAL.sam";
print NORMAL $head;
open ABNORMAL,">./$bamFile.ABNORMAL.sam";
print ABNORMAL $head;
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to sure that samtools is available \n";
# open BAM,"$bamFile" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to sure that samtools is available \n";

my %hash;
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
    #统计多处比对
    if($getSecondaryFlag || $getSupplementaryAlignmentFlag){
        ++ $count_multi;
        print MULTI "$_\n";
        next;
    }
    # next if($getSecondaryFlag || $getSupplementaryAlignmentFlag);
    ++ $totalReads;
    #统计带XA标签的
    if(/XA/){
        ++ $count_xa;
        print XA "$_\n";
        next;
    }
    #统计UNmap数
    if($getReadUnmappedFlag || $getMateUnmappedFlag){
        ++ $count_unmap;
        #这里算的是reads数 统计insert 时记得除以二
        print UMAP "$_\n";
        next;
    }
    # next if ($getReadUnmappedFlag || $getMateUnmappedFlag);

    # next if(!$getReadPairedFlag || $binFlag[8] || $getSupplementaryAlignmentFlag || $getReadFailsVendorQualityCheckFlag || $getReadUnmappedFlag);
        # 存入reads
        if($getFirstPairFlag){
            $hash{$line[0]}{data}{r1}=$_;
        }else{
            $hash{$line[0]}{data}{r2}=$_;
        }
        #判断方向
        my $distance;
        if($binFlag[4] == $binFlag[5]){
            $distance = 0;
        }else{
            my $positiveStart;
            my $negativeStart;
            my $matchLength;
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

        #判断是否正常
        if(!exists $hash{$line[0]}{normal}){
            if($cigar!~/S/ && abs $insertSize <= $max && $distance==1){
                $hash{$line[0]}{normal} = 1;
            }else{
                $hash{$line[0]}{normal} = 2
            }
        }elsif ($hash{$line[0]}{normal}==1){
            if($cigar!~/S/ && abs $insertSize <= $max && $distance==1){
                ++ $count_normal;
                print NORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
            }else{
                ++ $count_abnormal;
                print ABNORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";

            }
        }else{
            ++ $count_abnormal;
            print ABNORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
        }
}
my $total_insert = $totalReads/2;
my $pct_normal = sprintf "%.2f%%",100 * $count_normal/$total_insert;
my $pct_abnormal = sprintf "%.2f%%",100 * $count_abnormal/$total_insert;
my $pct_unmap = sprintf "%.2f%%",100 * $count_unmap/$total_insert;
my $pct_multi = sprintf "%.2f%%",100 * $count_multi/$total_insert;
print COUNT "COUNT_ALL_READS\tCOUNT_ALL_INSERT\tNORMAL_INSERT\tABNORMAL_INSERT\tUNMAP_INSERT\tSECONDERY\tXATag\tPCT_NORMAL_INSERT\tPCT_ABNORMAL_INSERT\tPCT_UNMAP_INSERT\n";
print COUNT "$totalReads\t$total_insert\t$count_normal\t$count_abnormal\t$count_unmap\t$count_multi\t$count_xa\t$pct_normal\t$pct_abnormal\t$pct_unmap\n";