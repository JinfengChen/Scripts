my @word=(1,2,3,4);

print "@word\n";

my $ref=\@word;

for(my $i=0;$i<@$ref;$i++){
   print "$i\n";
}
