=header
perl $0 test.fa > test.sam.header

@HD     VN:1.4  SO:coordinate
@SQ     SN:chr01        LN:38770143
@SQ     SN:chr02        LN:32956989
@SQ     SN:chr03        LN:34146764
@SQ     SN:chr04        LN:28064511
@SQ     SN:chr05        LN:24671766
@SQ     SN:chr06        LN:26531283
@SQ     SN:chr07        LN:25188480
@SQ     SN:chr08        LN:23957623
@SQ     SN:chr09        LN:19575785
@SQ     SN:chr10        LN:19599887
@SQ     SN:chr11        LN:24091154
@SQ     SN:chr12        LN:21977684

=cut

fasta2sam($ARGV[0]);

sub fasta2sam
{
$/=">";
my %hash;
my ($file)=@_;
print "\@HD	VN:1.4	SO:coordinate\n";
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $len = length $seq;
    print "\@SQ	SN:$head\tLN:$len\n";
}
close IN;
$/="\n";
return \%hash;
}

