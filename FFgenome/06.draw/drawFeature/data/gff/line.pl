#!/usr/bin/perl
my $gff=$ARGV[0];
chomp $gff;
my $readn;
my @all=<$gff/chr*>;
foreach (@all){
   my @read=split(" ",`wc -l $_`);
   $readn+=$read[0];
}
print "$readn\n";
