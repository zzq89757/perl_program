use strict;
my $pt = shift;
my $mt = shift;
my @pt;
my @mt;
open IN,$pt;
while(<IN>){
	chomp;
	next if (/^#/ || !/\S+/ || /^C/);
	my @line = split /\s+/;
	my $pct = sprintf "%.6f",$line[20];
	push @pt,$pct;
	
}
close IN;
open IN,$mt;
while(<IN>){
	chomp;
	next if (!/^P/);
	my @line = split/\s+/;
	push @mt,$line[1];
	
}
print "$mt checked pass!\n" if(@mt ~~ @pt);
