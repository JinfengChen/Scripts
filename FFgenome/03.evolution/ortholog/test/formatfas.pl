my %hash;
$/=">";
while(glob("*.fas")){
    print "$_\n";
    my $file=$_;
    my $spec;
    if ($_=~/(.*)\.fas/){
       $spec=$1;
    }
    open IN, "$file" or die "infile erros";
    open OUT, ">$spec\.named" or die "out erros";
    open OUT1, ">$spec\.length" or die "out2 erros";
    while(<IN>){
        my @unit=split("\n",$_);
        my $head=shift @unit;
        my $seq=join("",@unit);
        $seq=~s/>//;
        my @word=split(" ", $head);
        my $name=$word[0];
        
        my $length=length $seq;
       if ($name=~/\w+/){
        unless(exists $hash{$name}){
        $hash{$name}=$length;
        }else{
        print "$name\n";
        }
  
        print OUT ">$name\n$seq\n";
        print OUT1 "$name\t$length\n";
        }    
    }
    close IN;
    close OUT;
    close OUT1; 
}
