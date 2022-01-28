my $bam = shift;
my $r1 = shift;
my $r2 = shift;
my %hash;
open IN,"samtools view $bam|";
open R1,">$r1";
open R2,">$r2";
while(<IN>){
    chomp;
    my @line = split /\t/;
    my $name = $line[0];
    my $flag = $line[1];
    my $chr = $line[2];
    my $position = $line[3];
    my $cigar = $line[5];
    my $matePosition = $line[7];
    my $seq = $line[9];
    my $binnaryFlagString = sprintf("%b",$flag);
    my @binFlag =split //,$binnaryFlagString;
    @binFlag = reverse @binFlag;
    my $getFirstPairFlag = $binFlag[6];
  if(!exists $hash{count}{$name}){
    if($cigar=~/S/){
      if($getFirstPairFlag){
        $hash{data}{$name}{R1} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n";
        $hash{count}{$name} = 1;  
      }else{
        $hash{data}{$name}{R2} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n";
        $hash{count}{$name} = 1;  
      }
    }else{
      if($getFirstPairFlag){
        $hash{data}{$name}{R1} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n"; 
      }else{
        $hash{data}{$name}{R2} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n";
      }
        $hash{count}{$name} = 0;
    }
  }elsif($hash{count}{$name}==1){
      if($getFirstPairFlag){
        $hash{data}{$name}{R1} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n"; 
      }else{
        $hash{data}{$name}{R2} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n"; 
      }
    delete $hash{count}{$name};
  }else{
  if($cigar!~/S/){
    delete $hash{data}{$name}{R1};
    delete $hash{data}{$name}{R2};
    delete $hash{count}{$name};
    next;
  }else{
      if($getFirstPairFlag){
        $hash{data}{$name}{R1} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n"; 
      }else{
        $hash{data}{$name}{R2} = '>'.$name.'_'.$chr.'_'.$position.'_'.$cigar."\n".$seq."\n"; 
      }
      delete $hash{count}{$name};
  }
  }
}
for my $key (keys %{$hash{data}}){
  print R1 $hash{data}{$key}{R1};
  print R2 $hash{data}{$key}{R2};
}