#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"orth1:s","orth2:s","help");


my $help=<<USAGE;
Find lost ortholog in two ortholog tables. Check the output gene in gene_family of ortholog_paralog_pipeline. 
This will help us to find out what's wrong with programs.
perl $0 --orth1 ./distance_BGI/IRGSP2OBRACH.distance.txt --orth2 ./distance_159_noPOPTR/OS2OB.distance.tx
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $ref1=readpair($opt{orth1});
my $ref2=readpair($opt{orth2});

my $lost1=findlost($ref1,$ref2,"Common1.txt");
my $lost2=findlost($ref2,$ref1,"Common2.txt");

writelost($lost1,$opt{orth2},"LostIn2.txt");
writelost($lost2,$opt{orth1},"LostIn1.txt");

#######
sub readpair
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    push (@{$hash{$unit[0]}}, "$unit[1]\t$unit[3]\t$unit[4]");
}
close IN;
return \%hash;
}

sub writelost
{
my ($lost,$ort,$file)=@_;
open OUT, ">$file" or die "$!";
print OUT "Lost ortholog in file $ort\n";
    foreach(keys %$lost){
       my $list=join("\t",@{$lost->{$_}});
       print OUT "$_\t$list\n";
    }
close OUT;
}



sub findlost
{
my ($ref1,$ref2,$file)=@_;
my %lost;
open OUT, ">$file" or die "$!";
foreach(keys %$ref1){
   unless (exists $ref2->{$_}){
      #print "$_\t$ref1->{$_}\n";
      $lost{$_}=$ref1->{$_};
   }else{
      my @temp1=sort @{$ref1->{$_}};
      my @temp2=sort @{$ref2->{$_}};
      my $list1=join("\t",@temp1);
      my $list2=join("\t",@temp2);
      my @unit1=split("\t",$temp1[0]);
      my @unit2=split("\t",$temp2[0]);
      next if ($unit1[0] eq $unit2[0]);
      print OUT "$_\t$list1\n";
      print OUT "$_\t$list2\n";
   }
}
close OUT;
return \%lost;
}
