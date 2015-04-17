#!/share/raid12/chenjinfeng/tools/perl/bin/perl
use Getopt::Long;
use Data::Dumper;
use Graph;
GetOptions(\%opt,"querygff:s","targetgff:s","blast:s","help");

my $help=<<USAGE;
An undirected graph with genes as nodes and protein similarites as edge weights was constructed.
Protein similarities were derived from pair-wise local Smith-Waterman alignments (blastp). An e-value <=1e-15 and a minimal alignment coverage of >= 70% of both protein sizes were required.

OB.gff and OS.gff are annotation for a DNA fregment not whole genome. So the result will be gene pair for a fregment.
If whole genome annotation are given, it is all right. 

Run: perl undirected.pl -querygff OB.gff -targetgff OS.gff -blast allvsall.blasttable

-querygff : merge all the region into one gff file
-targetgff: merge all the region into one gff file
-blast : table format result of blast, convert from BGI blastparser.pl.

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $graph=Graph::Undirected->new;
#my $genepair=blastpair($opt{blast});

#foreach (keys %$genepair){
#    my ($gene1,$gene2)=split("\t",$_);
#    my $identity=$genepair->{$_};
#    $graph->add_weighted_edge($gene1,$gene2,$identity);
#}
$graph=blastgraph($opt{blast},$graph);

my $queryID;
if ($opt{querygff}=~/OS/i){
   $queryID =TigrGFFid($opt{querygff});
}else{
   $queryID =GFFid($opt{querygff});
}
my $targetID;
if ($opt{targetgff}=~/OS/i){
   $targetID=TigrGFFid($opt{targetgff});
}else{
   $targetID=GFFid($opt{targetgff}); 
}
open OUT, ">genepair.txt" or die "$!";
for(my $i=0;$i<@$queryID;$i++){
    for(my $j=0;$j<@$targetID;$j++){
        if ($graph->has_edge($$queryID[$i],$$targetID[$j])){
             my $w=$graph->get_edge_weight($$queryID[$i],$$targetID[$j]);
             print OUT "$$queryID[$i]\t$$targetID[$j]\t$w\n";
        }
    }
}
close OUT;
############################sub function#########################################

sub GFFid
{
my ($gff)=@_;
my @id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2] eq "mRNA"){
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            push (@id,$1);
        }
    }

}
close IN;
return \@id;
}


sub TigrGFFid
{
my ($gff)=@_;
my @id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2] eq "mRNA"){
        if ($unit[8]=~/Alias=(LOC_Os\w+\.\d+)/){
            push (@id,$1);
        }
    }

}
close IN;
return \@id;
}



sub blastgraph
{
####read a blasttable file, return a ref of hash contain pair of homologous gene that have e-value <=1e-15,
####and coverage of > 70% for both protein sizes;
	my ($blast,$graph)=@_;
	my %record;
	open IN, "$blast" or die "$!";
	<IN>;
	while(<IN>){
		my @unit=split("\t",$_);
		if ($unit[13] > 1e-15){
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
        if ($qlencoverage >= 0.5 and $hlencoverage >= 0.5){
	#if ($qhspcoverage >= 0.3 and $hhspcoverage >= 0.3 and $qlencoverage >= 0.5 and $hlencoverage >= 0.5){
		#my $temp="$queryid\t$hitid";
		#$homologous{$temp}=$identity;
                #print "$queryid\t$hitid\t$identity\n";
                $graph->add_weighted_edge($queryid,$hitid,$identity);
	}
}
#return \%homologous;
return $graph;
}

