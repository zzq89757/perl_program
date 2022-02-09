my $raw = shift;
open IN,"$raw";
while(my $id = <IN>){
    my $seq = <IN>;
    <IN>;
    <IN>;
    print $id.$seq;
}