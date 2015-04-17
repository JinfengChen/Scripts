#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;
GetOptions (\%opt,"rpkm:s","body:s","project:s","help");


my $help=<<USAGE;
Group gene into quintiles by expression level of RPKM (cutoff=1).
Plot RPKM with methylation level.
perl $0 -rpkm FF.rpkm.tophat -body FF.gene.CG.level.part -project FF
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $reflevel=partlevel($opt{body});
my $refrpkm=rpkmexpr($opt{rpkm});
my @gene=sort {$refrpkm->{$a} <=> $refrpkm->{$b}} keys %$refrpkm;
my $number=int (@gene/50);
my $counter=0;
my %group;
my %genenum;
open OUT, ">$opt{project}.mCvsRPKM.sum" or die "$!";
foreach (@gene){
   #print "$_\t$refrpkm->{$_}\n";
   $counter++;
   my $index=int ($counter/$number);
   if (exists $group{$index}){
       if (exists $reflevel->{$_}){
          my $temp=$reflevel->{$_};
          my $ref=$group{$index};
          my $t1=$ref->[0]+$temp->[0];
          my $t2=$ref->[1]+$temp->[1];
          my $t3=$ref->[2]+$temp->[2];
          $group{$index}=[$t1,$t2,$t3];
          #print "$group{$index}->[0]\n";
          $genenum{$index}++;
       }
   }else{
       if (exists $reflevel->{$_}){
          my $temp=$reflevel->{$_}; 
          $group{$index}=[$temp->[0],$temp->[1],$temp->[2]];
          #print "$group{$index}->[0]\n";
          $genenum{$index}++;
       }
   }
}
my @temp=sort {$a <=> $b} keys %group;
foreach (@temp){
   my $ref=$group{$_};
   #print $ref;
   #print "$_\t$ref->[0]\t$ref->[1]\t$ref->[2]\n";
   my $up=$ref->[0]/$genenum{$_};
   my $body=$ref->[1]/$genenum{$_};
   my $down=$ref->[2]/$genenum{$_};
   print OUT "$up\t$body\t$down\n";
}
close OUT;
& draw();

#############
sub draw
{
open DRAW, ">$opt{project}.mCvsRPKM.r" or die "$!";
print DRAW<<"END.";
pdf("$opt{project}.mCvsRPKM.pdf")
matrix(scan("$opt{project}.mCvsRPKM.sum"),ncol=3,byrow=T)->x
par(las=1)
plot(x[,1],type="l",col="red",font.main=3,lwd=2,ylim=c(0,0.3),xlab="Expression",ylab="Methylation level")
points(x[,2],type="l",col="blue",lwd=2)
points(x[,3],type="l",col="lightblue",lwd=2)
dev.off()
END.
close DRAW;
`R --vanilla -q <$opt{project}.mCvsRPKM.r`;
}


sub partlevel
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    next if ($unit[1] == 0 and $unit[2] ==0);
    next if ($unit[3] == 0 and $unit[4] ==0);
    next if ($unit[5] == 0 and $unit[6] ==0);
    my $uplevel=$unit[2]/($unit[1]+$unit[2]);
    my $bodylevel=$unit[4]/($unit[3]+$unit[4]);
    my $downlevel=$unit[6]/($unit[5]+$unit[6]);
    $hash{$unit[0]}=[$uplevel,$bodylevel,$downlevel];
}
close IN;
return \%hash;
}

sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    unless ($unit[1]=~/\-/ or $unit[1] == 0){
       $hash{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}
 
