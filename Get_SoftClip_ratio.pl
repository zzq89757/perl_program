#!/usr/bin/perl -w
use strict;

die "perl $0 <bam file> <out file>\n" unless (@ARGV==2);

my $in=shift;
my $out=shift;

open (IN,"samtools view $in | ") or die $!;
open (OUT,">$out");

my %hash;
while(<IN>){
	chomp;
	my @a=split /\s+/;
	$hash{all}{$a[0]}=1;
	if($a[5]=~/S/ && $a[5]!~/H/){
		if(!exists $hash{SoftClip}{$a[0]}){
			$hash{SoftClip}{$a[0]}=1;
		}else{
			$hash{SoftClip}{$a[0]."_2"}=1;
		}
	}
}
close IN;

my $All=scalar keys %{$hash{all}};
my $SoftClip=scalar keys %{$hash{SoftClip}};

my $per=sprintf ("%.6f",(scalar $SoftClip/(2*scalar $All)));
print OUT "All\t$All\nnum_soft\t$SoftClip\nSoftClip\t$per\n";


