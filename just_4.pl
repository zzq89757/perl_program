my $str = "3S87M9S";
# my $s = $1 if($str=~/^(\d+)S|(\d+)S$/);
my @array = $str=~/(\d+)S/g;
if((@array == 1 && $array[0]>4) || (@array ==2 && ($array[0]>4 || $array[1]>4))){
  $hash{S4}{$line[0]}.=""
}
for my $item (@array){
  if($item>3){
    return 1;
  }else{

  }
}
if($array[0]>3||$array[1]>3){
  print "ass"
}
# print $s;
