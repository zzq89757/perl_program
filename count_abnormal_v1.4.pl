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
#--------------------------- about v1.2 -----------------------------------
#0.关于正常插入片段评价标准
#本版本采用picard默认使用的 median + 10 * MAD 
#1.直接提取picard信息内容
#2.由于多次使用picard算法对方向进行判断，该版本将其封装为函数，结构更加清晰
use strict;
use Getopt::Long;
use POSIX;
use Statistics::Descriptive;
use List::Util qw/sum/;
use feature qw(signatures);
no warnings qw(experimental::signatures);
sub usage {
        print STDERR << "USAGE";
Description     Count_abnormal By zzq
Author                  zhengzhiqiang\@vazyme.com
Data                    2022.1
Version                 V1.4
Usage                   perl $0 bamFile countFile


 e.g.:
         perl $0  /xxx/xxx.bam xxx.count
USAGE

exit 0
}
usage() if(@ARGV<2);
# picard的方向判断算法 除FR以外皆划入异常
sub getOrientation($getReadReverseFlag,$getMateReverseFlag,$cigar,$position,$matePosition,$insertSize){
  my $orientation;
  # tandem
  if($getReadReverseFlag == $getMateReverseFlag){
    $orientation = 0;
  }else{
    my $positiveStart;
    my $negativeStart;
    if($getReadReverseFlag){
        my @match = $cigar =~/(\d+)[DMNX]/g;
        my $matchLength = sum @match;
        $positiveStart = $matePosition;
        $negativeStart = $position + $matchLength -1;
    }else{
      $positiveStart = $position;
      $negativeStart = $position + $insertSize;
    }
    $orientation = $negativeStart > $positiveStart ? 1 : 0;
  }
  return $orientation;
}
my $bamFile = shift;
my $countFile = shift;
my $head = `samtools view -H $bamFile`;
my $count_xa = 0;
my $count_multi = 0;
my $count_abnormal = 0;
my $count_normal = 0;
my @insert;
my %hash;
`mkdir -m 777 out_v1.2`;
open COUNT,">out_v1.2/$countFile";
open XA,">out_v1.2/$bamFile.XA.sam";
print XA $head;
open MULTI,">out_v1.2/$bamFile.MULTI.sam";
print MULTI $head;
open UMAP,">out_v1.2/$bamFile.UMAP.sam";
print UMAP $head;
open NORMAL,">out_v1.2/$bamFile.NORMAL.sam";
print NORMAL $head;
open ABNORMAL,">out_v1.2/$bamFile.ABNORMAL.sam";
print ABNORMAL $head;
print strftime("get max insert_size: %H:%M:%S\n", localtime);
# Insert Size #
my ($median,$mad,$mean);
if ( -e "/mnt/data/Users/zzq/Abnormal_align/zzq/cut/g4blast/test/pi_out/P6-1.CollectInsertSizeMetrics"){
	open (IN,"/mnt/data/Users/zzq/Abnormal_align/zzq/cut/g4blast/test/pi_out/P6-1.CollectInsertSizeMetrics");
	while(<IN>){
		chomp;
        next if(/^#/);
		if($_=~/^MEDIAN_INSERT_SIZE/){
			$_=<IN>;chomp;
			my @data=split("\t",$_);
            $median = $data[0];
            $mad = $data[1];
            $mean = $data[4];
            last;
		}
	}
	close IN;
}
my $max_insertSize = $median + 10 * $mad;
print "max insertSize is $max_insertSize\n";
# ESTIMATED LIBRARY SIZE #
my ($librarySize);
if ( -e "/mnt/data/Users/zzq/Abnormal_align/zzq/cut/g4blast/test/pi_out/P6-1.EstimateLibraryComplexity"){
	open (IN,"/mnt/data/Users/zzq/Abnormal_align/zzq/cut/g4blast/test/pi_out/P6-1.EstimateLibraryComplexity");
	while(<IN>){
		if($_=~/^LIBRARY\s+/){
			$_=<IN>;chomp;
			my @data=split(/\t/,$_);
            $librarySize = $data[9];
			# print OUT "Est_Library_Size\t$data[9]\n";
#			print OUT "ESTIMATED LIBRARY SIZE\t$data[9]\n";
		}
	}
	close IN;
}
 # 调用Samtools stats 进行各项指标统计
my ($clean_reads,$mapped_reads,$dup_reads,$clean_bases,$mapped_bases,$mismatch_bases,$uniq_reads,$mean_insert_size);
open IN,"samtools stats $bamFile|grep ^SN | cut -f 2-|";
while(<IN>){
  chomp;
  		if (/raw total sequences:\s+([0-9]+)/){
			$clean_reads = $1;
		}elsif (/reads mapped:\s+([0-9]+)/){
			$mapped_reads = $1;
		}elsif (/reads duplicated:\s+([0-9]+)/){
			$dup_reads = $1;
		}elsif (/total length:\s+([0-9]+)/){
			$clean_bases = $1;
		}elsif (/bases mapped \(cigar\):\s+([0-9]+)/){
			$mapped_bases = $1;
		}elsif (/mismatches:\s+([0-9]+)/){
			$mismatch_bases = $1;
		}elsif (/insert size average:\s+(\S+)/){
            $mean_insert_size = $1;
        }
}

print strftime(" read bamFile begin: %H:%M:%S\n", localtime);
open BAM,"samtools view $bamFile|" or  die "ERROR:Can't open file with samtools,please check your file and run 'samtools' to make sure that samtools is available \n";
# 读bam 得到统计表并拆分bam

while(<BAM>){
    chomp;
    my @line = split /\t/;
    my $flag = $line[1];
    my $position = $line[3];
    my $mappingQuality = $line[4];
    my $cigar = $line[5];

    my $matePosition = $line[7];
    my $insertSize = $line[8];
    my $seq = $line[9];
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
    my $getDupFlag = $binFlag[10];
    my $getSupplementaryAlignmentFlag = $binFlag[11];
    my $orientation;
    # 统计多处比对
    if($getSecondaryFlag || $getSupplementaryAlignmentFlag){
        ++ $count_multi;
        print MULTI "$_\n";
        next;
    }
    # 统计UNmap insert数
    if($getReadUnmappedFlag || $getMateUnmappedFlag){
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
        $orientation = getOrientation($getReadReverseFlag,$getMateReverseFlag,$cigar,$position,$matePosition,$insertSize);

        #判断是否正常
        if(!exists $hash{$line[0]}{normal}){
            if($cigar!~/S/ && abs $insertSize <= $max_insertSize && $orientation == 1){
                $hash{$line[0]}{normal} = 1;
            }else{
                $hash{$line[0]}{normal} = 2
            }
        }elsif ($hash{$line[0]}{normal}==1){
            if($cigar!~/S/ && abs $insertSize <= $max_insertSize && $orientation == 1){
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
print strftime(" read bamFile complete: %H:%M:%S\n", localtime);


print strftime("Classify Abnormal begin: %H:%M:%S\n", localtime);
#对异常情况进行分类
my $count_sa;
my $count_insert;
my $count_contig;
my $count_notfr;
my $count_others;
my $count_chimeras;
open INSERT,">out_v1.2/$bamFile.INSERT.sam";
print INSERT $head;
open CONTIG,">out_v1.2/$bamFile.CONTIG.sam";
print CONTIG $head;
open SA,">out_v1.2/$bamFile.SA.sam";
print SA $head;
open FR,">out_v1.2/$bamFile.FR.sam";
print FR $head;
open OTHERS,">out_v1.2/$bamFile.OTHERS.sam" or die ">out_v1.2/$bamFile.OTHERS.sam";
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
        delete $hash{abnormal}{$line1[0]};
    }else{
        ++ $count_chimeras;
        delete $hash{abnormal}{$line1[0]};
    }
}
print strftime("Classify Abnormal end: %H:%M:%S\n", localtime);
my $total_insert = $clean_reads/2;
my $pct_normal = sprintf "%.2f%%",100 * $count_normal/$total_insert;
my $pct_abnormal = sprintf "%.2f%%",100 * $count_abnormal/$total_insert;
my $pct_multi = sprintf "%.2f%%",100 * $count_multi/$total_insert;
my $pct_chimeras=sprintf "%.2f%%",100 * $count_chimeras/$clean_reads;
my $mapping_rate=sprintf("%.2f",$mapped_reads/$clean_reads*100);
my $dup_rate=sprintf("%.2f",$dup_reads/$mapped_reads*100);
my $uniq_rate=sprintf("%.2f",$uniq_reads/$mapped_reads*100);
my $mismatch_rate=sprintf("%.2f",$mismatch_bases/$mapped_bases*100);
print COUNT "CLEAN_READS\t$clean_reads\nCLEAN_BASES\t$clean_bases\nMAPPED_READS\t$mapped_reads\nMAPPED_BASES\t$mapped_bases\nMAPPED_RATE\t$mapping_rate\nDUP_READS\t$dup_reads\nDUP_RATE\t$dup_rate\nMISMATCH_BASES\t$mismatch_bases\nMISMATCH_RATE\t$mismatch_rate\nMEDIAN_INSERT_SIZE\t$median\nMEAN_INSERT_SIZE\t$mean\nLIBRARY_Complexity\t$librarySize\n\nCOUNT_ALL_INSERT\t$total_insert\nCOUNT_NORMAL_INSERT\t$count_normal($pct_normal)\nCOUNT_ABNORMAL_INSERT\t$count_abnormal($pct_abnormal)\nCOUNT_SECONDERY\t$count_multi\nCOUNT_XATag_INSERT\t$count_xa\nCOUNT_CHIMERAS_INSERT\t$count_chimeras($pct_chimeras)\n\nCOUNT_SATag_INSERT\t$count_sa\nCOUNT_NOT_FR\t$count_notfr\nCOUNT_OVERSIZE(over than $max_insertSize)\t$count_insert\nCOUNT_DIFF_CONTIG\t$count_contig\nCOUNT_OTHERS\t$count_others";