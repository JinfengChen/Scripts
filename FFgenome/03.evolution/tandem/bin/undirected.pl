#!/share/raid12/chenjinfeng/tools/perl/bin/perl
use Getopt::Long;
use Data::Dumper;
use Boost::Graph;
use Graph;
GetOptions(\%opt,"gff:s","blast:s","help");

my $help=<<USAGE;
An undirected graph with genes as nodes and protein similarites as edge weights was constructed.
Protein similarities were derived from pair-wise local Smith-Waterman alignments (blastp). An e-value <=1e-15 and a minimal alignment coverage of >= 70% of both protein sizes were required. Edges connecing genes that were more than 9 genes distant from each other in the genome were removed and tandem clusters were retrieved as connected groups from the resulting graph.

Run: perl undirected.pl -gff gene.gff -blast allvsall.blasttable

-gff   : gff annotation file
-blast : table format result of blast, convert from BGI blastparser.pl.

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

#my $graph =new Boost::Graph (directed=>0, net_name=>'Tandem',net_id=>1000);
#$graph->add_edge(node1=>'gene1',node2=>'gene2',weight=>1.2,edge=>'cluster');
#my $node=$graph->get_nodes();
#print Dumper($node);

my $graph=Graph::Undirected->new;
my $genepair=blastpair($opt{blast});
our %gene2pos;
our %gene2chr;

foreach (keys %$genepair){
    #print "$_\t$genepair->{$_}","\n";
    my ($gene1,$gene2)=split("\t",$_);
    my $identity=$genepair->{$_};
    $graph->add_weighted_edge($gene1,$gene2,$identity);
}
#$genepair="";
& gffinf($opt{gff});
my $nodenum=$graph->vertices(); 
my $edgenum=$graph->edges();
print "Node: $nodenum\nEdge: $edgenum\n";

######cut off weight(identity) and non tandem repeat for each edge#######
my @edges=$graph->edges();
foreach(@edges){
   #print "@{$_}\t";
   my $w=$graph->get_edge_weight($$_[0],$$_[1]);
   #print "$w\n";
   #if ($w < 0.3){
   #   $graph->delete_edge($$_[0],$$_[1]);
   #}
   unless (tandem($$_[0],$$_[1])){
      $graph->delete_edge($$_[0],$$_[1]);
   }
}
#%gene2chr={};
#%gene2pos={};
#########################################################################
=pod
###############delete single isolated node form graph####################
my @node=$graph->vertices();
foreach(@node){
   if ($graph->is_isolated_vertex($_)){
       $graph->delete_vertex($_);
       print "$_\n";
   }
}

#########################################################################
=cut

##############output cluster############################################
$nodenum=$graph->vertices();
$edgenum=$graph->edges();
print "Node: $nodenum\nEdge: $edgenum\n";
@cluster=$graph->connected_components();
my $tandem;
my $cluster;
open OUT, ">tandem_repeat_gene.txt";
for(my $i=0;$i<@cluster;$i++){
    if (@{$cluster[$i]} > 1){
       $cluster++;
       $tandem+=@{$cluster[$i]};
       print OUT "@{$cluster[$i]}\n";
    }
}
close OUT;
print "Tandem Cluster: $cluster\n";
print "Tandem Gene: $tandem\n";

###sub function############################
sub tandem
{
#### check if two gene are in tandem repeat. 
#### return ture if they are in the same chromosome or scaffold, and do not seperate by more than 9 genes.
my ($gene1,$gene2)=@_;
if (!exists $gene2chr{$gene1} or !exists $gene2chr{$gene2}){
    print "$gene1 or $gene2 not found in GFF\n";
}
if ($gene2chr{$gene1} eq $gene2chr{$gene2}){
    my $chr=$gene2chr{$gene1};
    my $refhash=$gene2pos{$chr};
    my $intergene=abs ($refhash->{$gene2}-$refhash->{$gene1});
    #if ($gene1 eq "Bradi3g28100.1" or $gene2 eq "Bradi3g28100.1"){
    #     print "$gene1\t$gene2\t$refhash->{$gene1}\t$refhash->{$gene2}\t$intergene\n";
    #}
    if ($intergene > 9){
       return 0;
    }else{
       return 1;
    } 
}else{
   return 0;
}
}

sub gffinf
{
###parse gff file store gene order information for each chr in hash
my ($gff)=@_;
#my %gene2pos;  ##hash to store gene position on scaffold or chromosome
#my %gene2chr;  ##hash to store chromosome that gene located on
my $seq;   ##Scaffold
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $gene2chr{$id}=$seq;
        if (exists $gene2pos{$seq}){
            my $refhash=$gene2pos{$seq};
            $refhash->{$id}=$unit[3];
            $gene2pos{$seq}=$refhash;
        }else{
            my %hash;
            $hash{$id}=$unit[3];
            $gene2pos{$seq}=\%hash;
        }
    }

}
close IN;

foreach(keys %gene2pos){
     my $seq=$_;
     my $refhash=$gene2pos{$_};
     my $counter;
     foreach(sort {$refhash->{$a} <=> $refhash->{$b}} keys %$refhash){
         $counter++;
         $refhash->{$_}=$counter;
         #print "$_\t$refhash->{$_}\n";          
     }
     $gene2pos{$seq}=$refhash;
}

#return (\%gene2pos,\%gene2chr);
}







sub blastpair
{
####read a blasttable file, return a ref of hash contain pair of homologous gene that have e-value <=1e-15,
####and coverage of > 70% for both protein sizes;
my ($blast)=@_;
my %record;
open IN, "$blast" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   if ($unit[13] > 1e-5){
       next;
   }
   my $pair=$unit[0]."-".$unit[4];
   if (exists $record{$pair}){
        my $refarray=$record{$pair};
	push (@$refarray,[@unit]);
	$record{$pair}=$refarray;
   }else{  
	my @array;
	push (@array,[@unit]);
        $record{$pair}=\@array;
   }
}
close IN;

my %homologous;
foreach (keys %record){
      #print "$_\n";
      #print Dumper($record{$_}),"\n";
      my $refarray=$record{$_};
      my $num=@$refarray;
      #print "$num\n";
      my $hsp=1;      
      my $queryid =$$refarray[0][0];
      my $hitid   =$$refarray[0][4];
      my $bitscore=$$refarray[0][12];
      my $identity=$$refarray[0][8];
      my $querylen=$$refarray[0][1];
      my $hitlen  =$$refarray[0][5];
      my $longestq=$$refarray[0][3]-$$refarray[0][2]+1;
      my $longesth=$$refarray[0][7]-$$refarray[0][6]+1;
      my $matchq  =$$refarray[0][3]-$$refarray[0][2]+1;
      my $matchh  =$$refarray[0][7]-$$refarray[0][6]+1;
      if ($num > 1){  
          my $hspstartq=$$refarray[0][2];
	  my $hspstarth=$$refarray[0][6];
	  my $hspendq  =$$refarray[0][3];
	  my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
	              next;
		}else{
                      $hsp++;
		      $hspendq=$$refarray[$i][3];
		      $hspendh=$$refarray[$i][7];
		      $bitscore+=$$refarray[$i][12];
                      $identity+=$$refarray[$i][8];
		      $matchq  +=$$refarray[$i][3]-$$refarray[$i][2]+1;
		      $matchh  +=$$refarray[$i][7]-$$refarray[$i][6]+1;
		} 
          } 
	  $longestq=$hspendq-$hspstartq+1;
	  $longesth=$hspendh-$hspstarth+1;
          $identity=$identity/$hsp;
      }
      my $qhspcoverage=$matchq/$querylen;
      my $hhspcoverage=$matchh/$hitlen;
      my $qlencoverage=$longestq/$querylen;
      my $hlencoverage=$longesth/$hitlen;
      #if ($qhspcoverage >= 0.7 and $hhspcoverage >= 0.7){
      #if ($qlencoverage >= 0.7 and $hlencoverage >= 0.7){
      if ($qhspcoverage >= 0 and $hhspcoverage >= 0 and $qlencoverage >= 0 and $hlencoverage >= 0){
          my $temp="$queryid\t$hitid";
          #print "$temp\n";
          $homologous{$temp}=$identity;
      }
}
return \%homologous;
}
