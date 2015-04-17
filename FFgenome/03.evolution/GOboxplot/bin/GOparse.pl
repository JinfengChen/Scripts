#!/usr/bin/perl
use Getopt::Long;
use GO::Parser;
use Data::Dumper;

GetOptions (\%opt,"obo:s","table:s","wego:s","column:s","level:s","project:s","help");


my $help=<<USAGE;
Parse gene_ontology.obo file, store the GO term of level 3 or 2 in hash.
For each GO term in hash, we calculate the number of gene acciocated with this term and statistics for values in column.
  
Example table.txt:
OBR_GLEAN_10001750      0.0855  0.5877  0.1454  45.20769231
OBR_GLEAN_10027965      0.1509  0.5122  0.2945  39.4

Result:
Level	GO term		Namespace		Annoation	Gene number   Mean    Median   Value             Gene
3       GO:0043025      cellular_component      cell body	123	      0.4     0.5      0.02,0.023,0.04   gene,gene

All: level 0
Molecular function, Biological process, Cellular component: level 1
......
If there are multi path to top we choose the shortest one.

Run: perl GOlevel.pl -obo ../input/gene_ontology.obo -wego OB.wego -table table.txt -column 2 -level 3 -project Ks
-obo: gene_ontology.obo
-wego: *.wego containing GO information
-table:  Gene information table
-column: 2 mean column 3 will be summarized to output
-level: GO term level to summarize
-project: prefix name for result file, such as OB2OS_Ks

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my ($refgene2go)=gene2go($opt{wego});
my $parser=new GO::Parser ({format=>'obo_text',
                            handler=>'obj'});

$parser->parse ("$opt{obo}");
my $graph=$parser->handler->graph;

#########store all level 3 term in hash###################
my %levelGO;
my $node=$graph->get_all_nodes();
foreach(@$node){
    my $term=$graph->get_term($_);
    my $id=$term->acc;
    my $name =$term->name;
    my $namespace=$term->namespace;
    #print "$id\t$name\n";
    my $path=$graph->paths_to_top($id);
    my @level;
    foreach(@$path){
      my $len=$_->length+1;
      push (@level,$len);
    }
    @level=sort {$a <=> $b} @level;
    my $level=shift @level;
    if ($level == $opt{level}){
      #print "$level\t$id\t$namespace\t$name\n";
      $levelGO{$id}=$level;
    }
}
###################################################


my %go2gene;
my %go2num;
my $gtotal;
my %gcount;
open IN, "$opt{table}" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_ eq "");
    $gtotal++;
    my @unit=split("\t",$_);
    if (exists $refgene2go->{$unit[0]}){
        #print "Gene:\t$unit[0]\n";
        my $refarray=$refgene2go->{$unit[0]};
        foreach(@$refarray){
                #print "GO:\t$_\n";
                my $path=$graph->paths_to_top($_);
                foreach(@$path){
                   my @temp;
                   #my $ref=$_->link_list();
                   my $ref=$_->term_list;
                   foreach(@$ref){
                       #next if ($_ =~/is_a/ or $_ =~/part_of/);
                       #print Dumper($_),"\n";
                       my $id=$_->acc;
                       my $name=$_->name;
                       #print "$id\t$name\n";
                       push (@temp,$id);
                   }
                   if (@temp >= $opt{level}){
                       $gcount{$unit[0]}=1;
                       my $index=@temp-$opt{level};
                       if (exists $go2gene{$temp[$index]}){  
                           my $refgene=$go2gene{$temp[$index]};
                           my $refnum =$go2num{$temp[$index]};
                           unless(exists $refgene->{$unit[0]}){  ## to keep one gene count only once for a GO term
                             $refgene->{$unit[0]}=1;
                             push (@$refnum,$unit[$opt{column}]);
                             $go2gene{$temp[$index]}=$refgene;
                             $go2num{$temp[$index]}=$refnum;
                             #print "Garray:@$refgene\nNarray:@$refnum\n";
                           }
                       }else{
                           my %gene;
                           my @num;
                           $gene{$unit[0]}=1;
                           push (@num,$unit[$opt{column}]);
                           $go2gene{$temp[$index]}=\%gene;
                           $go2num{$temp[$index]}=\@num;
                           #print "Garray:@gene\nNarray:@num\n";
                       }
                   }else{
                       print "@temp is up on Level $opt{level}\n";
                   }
                }
                #print "\n";
        }
    }else{
        print "GO annotation is not avaliable for term $unit[0]\n";
    }
}
close IN;
my $ngcount=keys %gcount;
print "Total Gene: $gtotal\n";
print "Count Gene: $ngcount\n";
###########################################################
open OUT, ">$opt{project}.GO.result" or die "$!";
foreach(sort keys %go2gene){
     my $ref1=$go2gene{$_};
     my $ref2=$go2num{$_};
     my $gnumber=keys %$ref1;
     my $glist=join(",",keys %$ref1);
     my $nlist=join(",",@$ref2);
     my ($mean,$median)=nstat(@$ref2);
     #print "$mean\t$median\n";
     my $term=$graph->get_term($_);
     my $id=$term->acc;
     my $name =$term->name;
     my $namespace=$term->namespace;
     my $level=$levelGO{$_};
     print OUT "$level\t$_\t$namespace\t$name\t$gnumber\t$mean\t$median\t$nlist\t$glist\n";
}
close OUT;

########sub function##############
sub gene2go{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
  chomp $_;
  next if ($_ eq "");
  my @unit=split("\t",$_);
  my $gene=shift @unit;
  $hash{$gene}=\@unit;
}
close IN;
return (\%hash);
}



sub nstat
{
my (@num)=@_;
@num=sort {$a<=>$b} @num;
my $loop=0;
my $total=0;
my $add_square=0;
foreach  (@num) {
        $total+=$_;
        $add_square+=$_*$_;
        $loop++;
}
my $number=@num;
my $min=$num[0];
my $max=$num[$number-1];
my $mean=$total/$number;
my $median=$num[int $number/2];
#$SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
return ($mean,$median);
}
