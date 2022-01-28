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
#--------------------------- about v1.1 -----------------------------------
#0.关于正常插入片段评价标准
#本版本采用picard默认使用的 median + 10 * MAD 
#1.单个脚本即可将MAD算出并同时统计出各种情况的数目并拆分bam文件但需要读取两遍bam
#2.尽管本程序读取了两次bam文件，但相较上一版本效率仍提升一倍多
#3.对程序整体速度进行优化，由原来的 20min/GB 提升至 10min/GB
#4.由于多次使用picard算法对方向进行判断，该版本将其封装为函数，结构更加清晰
use strict;
use Getopt::Long;
use POSIX;
use Statistics::Descriptive;
use feature qw(signatures);
no warnings qw(experimental::signatures);
sub usage {
        print STDERR << "USAGE";
Description     Count_abnormal By zzq
Author                  zhengzhiqiang\@vazyme.com
Data                    2021.12
Version                 V1.1
Usage                   perl $0 bamFile countFile


 e.g.:
         perl $0  /xxx/xxx.bam xxx.count
USAGE

exit 0
}
usage() if(@ARGV<2);
#picard的方向判断算法
sub getOrientation($getReadReverseFlag,$getMateReverseFlag,$cigar,$position,$matePosition,$insertSize){
  my $distance;
  if($getReadReverseFlag == $getMateReverseFlag){
    $distance = 0;
  }else{
    my $positiveStart;
    my $negativeStart;
    my $matchLength;
    if($getReadReverseFlag){
      my @match = $cigar =~/(\d+)[DMNX]/g;
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
  return $distance;
}
my $bamFile = shift;
my $countFile = shift;
my $head = `samtools view -H $bamFile`;
my $totalReads = 0;
my $count_xa = 0;
my $count_multi = 0;
my $count_unmap = 0;
my $count_abnormal = 0;
my $count_normal = 0;
my @insert;
my %hash;
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
print strftime("First time read bamFile begin: %H:%M:%S\n", localtime);
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to make sure that samtools is available \n";
#第一次读bam 得到绝对中位差
while(<BAM>){
    my @line = split /\t/;
    my $flag = $line[1];
    my $insertSize = abs $line[8];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getReadPairedFlag = $binFlag[0];
    my $getReadUnmappedFlag = $binFlag[2];
    my $getMateUnmappedFlag = $binFlag[3];
    my $getFirstPairFlag = $binFlag[6];
    my $getSecondaryFlag = $binFlag[8];
    my $getDupFlag = $binFlag[10];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    next if(!$getReadPairedFlag || $getSecondaryFlag || $getSupplementaryAlignmentFlag || $getDupFlag || $getReadUnmappedFlag || !$getFirstPairFlag || $getMateUnmappedFlag || $insertSize == 0);
    push (@insert,$insertSize);
}
close BAM;
print strftime("First time read bamFile complete，get max insertSize: %H:%M:%S\n", localtime);
my $record = Statistics::Descriptive::Full->new();
$record->add_data(\@insert);
my $median = $record->median();
my $mad = $record->median_absolute_deviation();
my $max_insertSize = $median + 10 * $mad;
print "max insertSize is $max_insertSize\n";
print strftime("Second time read bamFile begin: %H:%M:%S\n", localtime);
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to make sure that samtools is available \n";
#第二次读bam 得到统计表并拆分bam
while(<BAM>){
    chomp;
    my @line = split /\s+/;
    my $flag = $line[1];
    my $position = $line[3];
    my $mappingQuality = $line[4];
    my $cigar = $line[5];
    if($cigar=~/\*/){
        print "$line[0]\n";
    }
    my $matePosition = $line[7];
    my $insertSize = $line[8];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getReadPairedFlag = $binFlag[0];
    my $getReadUnmappedFlag = $binFlag[2];
    my $getMateUnmappedFlag = $binFlag[3];
    my $getReadReverseFlag = $binFlag[4];
    my $getMateReverseFlag = $binFlag[5];
    my $getFirstPairFlag = $binFlag[6];
    my $getSecondaryFlag = $binFlag[8];
    my $getReadFailsVendorQualityCheckFlag = $binFlag[9];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    my $distance;
    #统计多处比对
    if($getSecondaryFlag || $getSupplementaryAlignmentFlag){
        ++ $count_multi;
        print MULTI "$_\n";
        next;
    }
    ++ $totalReads;
    #统计UNmap数
    if($getReadUnmappedFlag || $getMateUnmappedFlag){
        ++ $count_unmap;
        print UMAP "$_\n";
        next;
    }
        # 存入reads
    if($getFirstPairFlag){
        $hash{$line[0]}{data}{r1} = $_;
    }else{
        $hash{$line[0]}{data}{r2} = $_;
    }
        #统计带XA标签的
    if(!exists $hash{$line[0]}{XA}){
        if(/XA:Z:/){
            $hash{$line[0]}{XA} = 1;
        }else{
            $hash{$line[0]}{XA} = 0
        }
    }elsif ($hash{$line[0]}{XA} == 1){
        ++ $count_xa;
        print XA $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
        delete $hash{$line[0]};
        next;
    }else{
        if(/XA:Z:/){
            ++ $count_xa;
            $hash{$line[0]}{XA} = 1;
            print XA $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
            delete $hash{$line[0]};
            next;
        }
    }
        #判断方向
        $distance = getOrientation($getReadReverseFlag,$getMateReverseFlag,$cigar,$position,$matePosition,$insertSize);

        #判断是否正常
        if(!exists $hash{$line[0]}{normal}){
            if($cigar!~/S/ && abs $insertSize <= $max_insertSize && $distance == 1){
                $hash{$line[0]}{normal} = 1;
            }else{
                $hash{$line[0]}{normal} = 2
            }
        }elsif ($hash{$line[0]}{normal}==1){
            if($cigar!~/S/ && abs $insertSize <= $max_insertSize && $distance == 1){
                ++ $count_normal;
                print NORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
                delete $hash{$line[0]};
            }else{
                ++ $count_abnormal;
                print ABNORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
            }
        }else{
            ++ $count_abnormal;
            print ABNORMAL $hash{$line[0]}{data}{r1}."\n".$hash{$line[0]}{data}{r2}."\n";
        }
}
print strftime("Second time read bamFile complete: %H:%M:%S\n", localtime);
my $total_insert = $totalReads/2;
my $pct_normal = sprintf "%.2f%%",100 * $count_normal/$total_insert;
my $pct_abnormal = sprintf "%.2f%%",100 * $count_abnormal/$total_insert;
my $pct_unmap = sprintf "%.2f%%",100 * $count_unmap/$total_insert;
my $pct_multi = sprintf "%.2f%%",100 * $count_multi/$total_insert;
my $count_unmap_insert = $count_unmap/2;

print strftime("Hash begin: %H:%M:%S\n", localtime);
#对异常情况进行分类
my $count_sa;
my $count_insert;
my $count_contig;
my $count_notfr;
my $count_others;
my $count_chimeras;
open INSERT,">./$bamFile.INSERT.sam";
print INSERT $head;
open CONTIG,">./$bamFile.CONTIG.sam";
print CONTIG $head;
open SA,">./$bamFile.SA.sam";
print SA $head;
open FR,">./$bamFile.FR.sam";
print FR $head;
open OTHERS,">./$bamFile.OTHERS.sam" or die "fuck\n";
print OTHERS $head;
foreach my $key (keys %hash){
    my $r1 = $hash{$key}{data}{r1};
    my $r2 = $hash{$key}{data}{r2};

    my @line1 = split /\t/,$r1;
    my $flag1 = $line1[1];
    my $position1 = $line1[3];
    my $cigar1 = $line1[5];
    my $matePosition1 = $line1[7];
    my $insertSize1 = $line1[8];
    my $binnaryFlagString1 = sprintf("%b",$flag1);
    my @binFlag1 = split //,$binnaryFlagString1;
    @binFlag1 = reverse @binFlag1;
    my $getReadReverseFlag1 = $binFlag1[4];
    my $getMateReverseFlag1 = $binFlag1[5];

    my @line2 = split /\t/,$r2;
    my $flag2 = $line2[1];
    my $position2 = $line2[3];
    my $cigar2 = $line2[5];
    my $matePosition2 = $line2[7];
    my $insertSize2 = $line2[8];
    my $binnaryFlagString2 = sprintf("%b",$flag2);
    my @binFlag2 = split //,$binnaryFlagString2;
    @binFlag2 = reverse @binFlag2;
    my $getReadReverseFlag2 = $binFlag2[4];
    my $getMateReverseFlag2 = $binFlag2[5];
    #SATag
    if($r1=~/SA:Z:/ || $r2=~/SA:Z:/){
        ++ $count_sa;
        $hash{abnormal}{$line1[0]} = 1;
        print SA "$r1\n$r2\n";
    }
    #判断方向 已封装函数 直接调用即可
    if(getOrientation($getReadReverseFlag1,$getMateReverseFlag1,$cigar1,$position1,$matePosition1,$insertSize1) == 0 || getOrientation($getReadReverseFlag2,$getMateReverseFlag2,$cigar2,$position2,$matePosition2,$insertSize2) == 0){
        ++$count_notfr;
        $hash{abnormal}{$line1[0]} = 1;
        print FR "$r1\n$r2\n";
    }
    #判断插入片段长度
    if(abs $line1[8] > $max_insertSize){
        ++ $count_insert;
        $hash{abnormal}{$line1[0]} = 1;
        print INSERT "$r1\n$r2\n";
    }
    #判断是否在同一染色体
    if($line1[6] !~ /=/ || $line2[6] !~ /=/){
        ++ $count_contig;
        $hash{abnormal}{$line1[0]} = 1;
        print CONTIG "$r1\n$r2\n";
    }
    #判断others
    if(!exists$hash{abnormal}{$line1[0]}){
        ++ $count_others;
        print OTHERS "$r1\n$r2\n";
    }else{
        ++ $count_chimeras;
    }
}
my $pct_chimeras=sprintf "%.2f%%",100 * $count_chimeras/$total_insert;
print strftime("Hash end: %H:%M:%S\n", localtime);
print COUNT "COUNT_ALL_INSERT\t$total_insert\nCOUNT_NORMAL_INSERT\t$count_normal($pct_normal)\nCOUNT_ABNORMAL_INSERT\t$count_abnormal($pct_abnormal)\nCOUNT_UNMAP_INSERT\t$count_unmap_insert($pct_unmap)\nCOUNT_SECONDERY\t$count_multi\nCOUNT_XATag_INSERT\t$count_xa\nCOUNT_CHIMERAS_INSERT\t$count_chimeras($pct_chimeras)\n\nCOUNT_SATag_INSERT\t$count_sa\nCOUNT_NOT_FR\t$count_notfr\nCOUNT_OVERSIZE(over than $max_insertSize)\t$count_insert\nCOUNT_DIFF_CONTIG\t$count_contig\nCOUNT_OTHERS\t$count_others";