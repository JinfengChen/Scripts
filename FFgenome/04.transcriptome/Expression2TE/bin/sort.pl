my $array=[0,5,2,6,7];
@$array=sort {$a <=> $b} @$array;
foreach(@$array){
   print "$_\n";
}
