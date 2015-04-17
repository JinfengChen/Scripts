#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"methyCG:s","methyCHG:s","methyCHH:s","cutoff:s","gff:s","project:s","help");


my $help=<<USAGE;
Split repeat GFF file into methylation/unmethylation by methylation level (uM if TE with methylation level < 10%).
--methyCG: REPEAT.CG.level.part file generate by methylevel.pl.
--methyCHG:
--methyCHH:
--gff: gff file of TE 
ID      UpC     UpmC    BodyC   BodymC  DownC   DownmC
TE063920        224     61      12      10      35      17
perl $0 --methyCG ../input/rice.REPEAT.CG.level.part --methyCHG ../input/rice.REPEAT.CHG.level.part --methyCHH ../input/rice.REPEAT.CHH.level.part --gff ../input/rice.repeat.gff --project rice
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{project} ||= "test";
$opt{cutoff} ||= 0.1;
my $refhashcg= TEmethylevel($opt{methyCG});
#& drawlevel($refhashcg,"CG");
& markTEbyMe($opt{gff},$refhashcg,"CG");

my $refhashchg= TEmethylevel($opt{methyCHG});
#& drawlevel($refhashchg,"CHG");
& markTEbyMe($opt{gff},$refhashchg,"CHG");

my $refhashchh= TEmethylevel($opt{methyCHH});
#& drawlevel($refhashchh,"CHH");
& markTEbyMe($opt{gff},$refhashchh,"CHH");

my %hash;
foreach(keys %$refhashcg){
     my $totalc=$refhashcg->{$_}->[0]+$refhashchg->{$_}->[0]+$refhashchh->{$_}->[0];
     my $totalmc=$refhashcg->{$_}->[1]+$refhashchg->{$_}->[1]+$refhashchh->{$_}->[1];
     $hash{$_}=[$totalc,$totalmc];
}
#& drawlevel(\%hash,"Methylation");
& markTEbyMe($opt{gff},\%hash,"Me");
& splitTEbyMe($opt{gff},\%hash);

#############################################
sub drawlevel{
my ($refhash,$title)=@_;
open OUT, ">$opt{project}.$title.4r" or die "$!";
foreach(keys %$refhash){
    my $level=$refhash->{$_}->[0] == 0 ? 0 : $refhash->{$_}->[1]/$refhash->{$_}->[0];
    print OUT "$_\t$level\n";
}
close OUT;

open OUT, ">$opt{project}.$title.r" or die "$!";
print OUT <<"END.";
read.table("$opt{project}.$title.4r") -> x
pdf("$opt{project}.$title.pdf")
hist(x[,2],breaks=20,col=2,xlim=c(0,1),xlab="Methylation Level",ylab="Frequency",main="$title")
dev.off()
END.
close OUT;

system ("cat $opt{project}.$title.r | R --vanilla --slave");
}



#########################################
sub TEmethylevel
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   next if ($_ =~/^$/);
   my @unit=split("\t",$_);
   #my $level=$unit[3] == 0 ? 0 : $unit[4]/$unit[3];
   #$hash{$unit[0]}=$level;
   $hash{$unit[0]}=[$unit[3],$unit[4]];
}
close IN;
return \%hash;
}

### Mark TE by methylation status, > 10% methylation.
sub markTEbyMe
{
my ($gff,$hash,$mark)=@_;
my $me;
if ($gff=~/(.*)\.gff/){
   $me=$1.".gff.$mark.status";
}
my $refgff=parseGFF($gff);
my ($total,$miss,$meth,$unmeth);
open OUT, ">$me" or die "$!";
foreach(keys %$refgff){
    $total++;
    if(exists $hash->{$_} and $hash->{$_}->[0] > 0){
        my $level=$hash->{$_}->[1]/$hash->{$_}->[0];
        if ($level >= $opt{cutoff}){
            $meth++;
            print OUT "$_\t$level\n";
        }else{
            $unmeth++;
            print OUT "$_\t$level\n";
        }
    }else{
        $miss++;
        print OUT "$_\tNA\n";
    }
}
close OUT;
print "Total TE Copy: $total\nNot Covered: $miss\nMethylation: $meth\nUnmethylation: $unmeth\n";
}



### split TE by methylation status, > 10% methylation.
sub splitTEbyMe
{
my ($gff,$hash)=@_;
my ($me,$unme);
if ($gff=~/(.*)\.gff/){
   $me=$1.".me.gff";
   $unme=$1.".unme.gff";
}
my $refgff=parseGFF($gff);
my ($total,$miss,$meth,$unmeth);
open OUT1, ">$me" or die "$!";
open OUT2, ">$unme" or die "$!";
foreach(keys %$refgff){
    $total++;
    if(exists $hash->{$_} and $hash->{$_}->[0] > 0){
        my $level=$hash->{$_}->[1]/$hash->{$_}->[0];
        if ($level >= $opt{cutoff}){
            $meth++;
            print OUT1 "$refgff->{$_}\n";
        }else{
            $unmeth++;
            print OUT2 "$refgff->{$_}\n";
        }
    }else{
        $miss++;
    }
}
close OUT2;
close OUT1;
print "Total TE Copy: $total\nNot Covered: $miss\nMethylation: $meth\nUnmethylation: $unmeth\n";
}

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/ or $_=~/^$/);
    $record=$_;
    if ($record=~/ID=(.*?);/){
            $id=$1;
    }
    #print "$id\t$record\n";
    $hash{$id}=$record;
}
close IN;
return \%hash;
}







