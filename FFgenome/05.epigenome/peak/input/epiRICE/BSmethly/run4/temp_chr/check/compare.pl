
my $hash1=readfile($ARGV[0]);
my $hash2=readfile($ARGV[1]);

foreach(keys %$hash1){
    if (exists $hash2->{$_}){
        #print "$_\n";
    }else{
        print "No\t$_\t$hash1->{$_}\n";
    }
}


sub readfile
{
my ($file)=@_;
my %hash;
open IN, "$file" or die;
while(<IN>){
    my @unit=split("\t",$_);
    my $index=$unit[0]."_".$unit[1];
    $hash{$index}=$_;
}
close IN;
return \%hash;
}
