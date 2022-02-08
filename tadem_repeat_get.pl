use strict;
# my $fa = shift;
# open IN,"$fa";
# while(my $id = <IN>){
#   my $seq = <IN>;
#   if()
# }
# my $str = "ACTAGAGAG";
my $str = "AGGGGG";
#匹配串联重复 但能同时匹配到poly结构 需要先过掉poly结构
if($str=~/([AGCT])\1{3}/){

}else{
  if($str=~/([ACGT]{2,3})\1{2}/){
  print "right";
}else{
  print "wrong";
}
}
