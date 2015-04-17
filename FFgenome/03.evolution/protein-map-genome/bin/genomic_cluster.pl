#!/usr/bin/perl

=head1 Usage
	perl $0 [option] gff_file/table_file
	--type  <str> set type of input file ,default='table'
	--overlap_length <int> set the cutoff for overlap length, recommend 100 bp.
	--overlap_percent <float> set the cutoff for overlap percent, recommend 0.1.

=cut

use strict;
use Getopt::Long;


my ($Overlap_length,$Overlap_percent);
my ($Help);
GetOptions(
	"overlap_length:i"=>\$Overlap_length,
	"overlap_percent:f"=>\$Overlap_percent,
	"help"=>\$Help,
);
die `pod2text $0` if(@ARGV==0 || $Help || !($Overlap_length || $Overlap_percent));
my $InFile=shift;
my %genes;

#table:
#At1g01010.1     Chr12   +       1138400 1146949  cds_len

Read_Table($InFile,\%genes);

foreach my $chr(sort keys %genes) {
	my @sort_genes= sort {$a->[3]<=>$b->[3]} @{$genes{$chr}};
	my @cluster;
	foreach my $gene (@sort_genes) {
		if (@cluster==0) {
			push @cluster,[@$gene];
		}else{
			my $overlap;
			if (defined $Overlap_length){
				$overlap=find_overlap(\@cluster,$gene,$Overlap_length,'length');
			}elsif(defined $Overlap_percent){
				$overlap=find_overlap(\@cluster,$gene,$Overlap_percent,'percent');
			}else{
				die "ERROR: input cutoff !\n";
			}	
			if ($overlap==1) {
				push @cluster,[@$gene];
			}else{
				my @noredundance=find_nonredundance(\@cluster);
				@cluster=();
				push @cluster,[@$gene];
				print join("\t",@noredundance)."\n";
			}
		}
	}
	my @noredundance=find_nonredundance(\@cluster);
	@cluster=();
	print join("\t",@noredundance)."\n";
}




#At1g01010.1     Chr12   +       1138400 1146949
sub Read_Table {
	my $InFile=shift;
	my $genes_p=shift;
	open IN,$InFile;
	while(<IN>) {
		chomp;
		my @c=split(/\t/,$_);
		@c[3,4]=($c[3]<$c[4])? @c[3,4]:@c[4,3];
		push @{$genes_p->{"$c[1]$c[2]"}},[@c];
	}
	close IN;	
}

sub find_overlap {
	my ($c_p,$g_p,$cutoff,$tag)=@_;
	my $gene_len=abs($g_p->[4]-$g_p->[3])+1;
	my $cutoff_len;
	my $overlap=0;
	foreach my $c_gene(@$c_p) {
		my $c_gene_start=$c_gene->[3];
		my $c_gene_end=$c_gene->[4];
		my $c_gene_len=$c_gene_end-$c_gene_start+1;
		if ($tag eq 'length'){
			$cutoff_len=$cutoff;
		}elsif($tag eq 'percent'){
			$cutoff_len=($gene_len<$c_gene_len)?$cutoff*$gene_len:$cutoff*$c_gene_len;
		}else{
			die "\nERROR: Unknow tag.\n";
		}
		if ($g_p->[4]<$c_gene_start||$g_p->[3]>$c_gene_end) {
			$overlap=0;
		}else{
			my $len=0;
			if($g_p->[4]<=$c_gene_end){
				$len=$gene_len;
			}else{
				$len=$c_gene_end-$g_p->[3]+1;
			}
			if ($len>=$cutoff_len){
				$overlap=1;
				last;
			}
		}
	}
	return $overlap;
}

sub find_nonredundance {
	my $cluster_p=shift;
	my $max_length=0;
	my @max_gene;
	my @redundance;
	my @noredundance;	
	foreach my $gene (@$cluster_p) {
		if ($gene->[5] < $max_length) {
			next;
		}else{
			$max_length=$gene->[5];
			@max_gene=@$gene;
		}
	}
	@noredundance=@max_gene;
	return @noredundance;
}
