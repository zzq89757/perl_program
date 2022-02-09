use strict;
my $bam = shift;
my ($clean_reads,$mapped_reads,$dup_reads,$clean_bases,$mapped_bases,$mismatch_bases,$uniq_reads);
open IN,"samtools stats $bam|grep ^SN | cut -f 2-|";
while(<IN>){
  chomp;
  		if (/raw total sequences:\s+([0-9]+)/){
			$clean_reads=$1;
		}elsif (/reads mapped:\s+([0-9]+)/){
			$mapped_reads=$1;
		}elsif (/reads duplicated:\s+([0-9]+)/){
			$dup_reads=$1;
		}elsif (/total length:\s+([0-9]+)/){
			$clean_bases=$1;
		}elsif (/bases mapped \(cigar\):\s+([0-9]+)/){
			$mapped_bases=$1;
		}elsif (/mismatches:\s+([0-9]+)/){
			$mismatch_bases=$1;
		}
}
	my $mapping_rate=sprintf("%.2f",$mapped_reads/$clean_reads*100);
	my $dup_rate=sprintf("%.2f",$dup_reads/$mapped_reads*100);
	my $uniq_rate=sprintf("%.2f",$uniq_reads/$mapped_reads*100);
	my $mismatch_rate=sprintf("%.2f",$mismatch_bases/$mapped_bases*100);
print "TOTAL_READS\t$clean_reads\nTOTAL_BASES\t$clean_bases\nMAPPED_READS\t$mapped_reads\nMAPPED_BASES\t$mapped_bases\nMAPPED_RATE\t$mapping_rate\nDUP_READS\t$dup_reads\nDUP_RATE\t$dup_rate\nMISMATCH_BASES\t$mismatch_bases\nMISMATCH_RATE\t$mismatch_rate\n";