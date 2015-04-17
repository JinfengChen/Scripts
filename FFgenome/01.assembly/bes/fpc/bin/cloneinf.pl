### read fpc file and get the information of each clone, output to a file
### usage: perl cloneinf.pl *.fpc

if (@ARGV<1){
      die "Usage: perl cloneinf.pl *.fpc \n";
      exit(0);    
}

my $file=$ARGV[0];
open INFILE, "$file" or die "can not open my file ";
while (<INFILE>){
    chomp $_;
    if ($_=~/\/\//){
       next;
    }elsif($_=~/BAC : \"(\w+)\"/){
       my $clone=$1;
       my $contig;
       my $leftp;
       my $rightp;
       my $left=<INFILE>;
       if ($left=~/Map \"(\w+)\" Ends Left (\d+)\.\d+/){
             $contig=$1;
             $leftp=$2;
       }else{
            next;
       }
       my $right=<INFILE>;
       if ($right=~/Map \"(\w+)\" Ends Right (\d+)\.\d+/){
             $rightp=$2;
       }else{
            next;
       }
       $contig=~/ctg(\d+)/;
       $contig=$1;
       print "$clone\t$contig\t$leftp\t$rightp\n";
    }
}
close INFILE;
