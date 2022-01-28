# use v5.20;
use feature qw(signatures);
no warnings qw(experimental::signatures);
sub test($a,$b,$c){
    # my($a,$b,$c)=shift;
     $b-$a+$c;
}
my $n = test(2,3,4);
print $n;