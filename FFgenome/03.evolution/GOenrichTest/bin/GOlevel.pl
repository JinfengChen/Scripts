#!/usr/bin/perl
use Getopt::Long;
use GO::Parser;
use Data::Dumper;

GetOptions (\%opt,"obo:s","blast2go:s","level:s","help");


my $help=<<USAGE;
Parse gene_ontology.obo file, then find level of each GO term in the blast2go GO enrichment test result file. 
Output record into file if it is at a given level.
All: level 0
Molecular function, Biological process, Cellular component: level 1
......
If there are multi path to top we choose the shortest one.

Run: perl GOlevel.pl -obo ../input/gene_ontology.obo -blast2go ../input/OB.GOenrichmentTest.txt -level 3
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my ($head,$list)=getlist($opt{blast2go});
my $parser=new GO::Parser ({format=>'obo_text',
                            handler=>'obj'});

$parser->parse ("$opt{obo}");
my $graph=$parser->handler->g;


print "$head\n";
foreach(keys %$list){
   my $go=$_;
   my $path=$graph->paths_to_top($go);
   my @level;
   foreach(@$path){
      my $len=$_->length+1;
      push (@level,$len);
   }
   @level=sort {$a <=> $b} @level;
   my $level=shift @level;
   if ($level == $opt{level}){
      print "$level\t$list->{$go}\n";
   }else{
      print "$level\t$list->{$go}\n";
   }
}


########sub function##############
sub getlist{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
my $head=<IN>;
chomp $head;
$head="Level\t".$head;
while(<IN>){
  chomp $_;
  next if ($_ eq "");
  my @unit=split("\t",$_);
  $hash{$unit[0]}=$_;
}
close IN;
return ($head,\%hash);
}

