#!/usr/bin/perl

=head1 Name

taxonomy_evolution.pl  --  draw phylogeny tree for a set of species

=head1 Description

This program is used to draw phlogeny tree of the given species, using NCBI taxomony database.
Note that it not restricted to species, but also can be genus or the other taxomony ranks.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage

  perl taxonomy_evolution.pl [options] <input_file>
  --outdir <str>   set the result directory
  --debug          output the middle result for checking error
  --verbose        output running progress information to screen  
  --help           output help information to screen  

=head1 Exmple

  perl ../bin/taxonomy_evolution.pl ../input/animals.list 
  perl ../bin/taxonomy_evolution.pl ../input/animals.list -verbose
  perl ../bin/taxonomy_evolution.pl ../input/animals.list -debug

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$Bin/../lib";
use File::Basename qw(basename dirname); 
use Data::Dumper;
use Tree::nhx_svg;


##get options from command line into variables and set default values
my ($Outdir,$Debug,$Verbose,$Help);
GetOptions(
	"outdir:s"=>\$Outdir,
	"debug:s"=>\$Debug,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Outdir =~ s/\/$//;
$Outdir ||= ".";
die `pod2text $0` if (@ARGV == 0 || $Help);

my $infile=shift;
my $file_core = basename($infile);

my %node;  #global
my %id_name; #global
my %id_parents; #global;  key: node id;  value:[1,1223,....,node id]
$id_name{0}="Empty_node"; #global
my $tree="No1de[&&NHX:T=1:S=root];"; #root node Tax_id is 1
my $empty="No0de[&&NHX:T=0:S=node_saver]";  #empty child node to save the parent node 

my %unfound; #store unfoud names or ids
my @born; #store nodes which have children 
my %name_disp;


print STDERR "\n\tInput file read\n" if(defined $Verbose);

&read_file($infile); #read the input file to costruct %id_name

print STDERR "\n\tNames or ids checked\n"  if(defined $Verbose);

&construct_id_parents();

print STDERR "\n\tId parental chain constructed\n"  if(defined $Verbose);

&construct_node_link();

&construct_node_anno();

print STDERR "\n\tNode structure constructed\n"  if(defined $Verbose) ;

&write_tree();

&draw_tree();

print STDERR "\n\tNHX tree written, SVG drawn\n"  if(defined $Verbose);

&disp_except() if(defined $Verbose);



####################################################
################# Main  Function ###################
####################################################





#read the input file and judge its types to construct %id_name
####################################################
sub read_file{
	my $infile=shift;
	open (IN, $infile) || die ("can not open $infile\n");
	my @lines;
	while (<IN>) {
		chomp;
		next unless($_);
		my ($species, $display_name) = (split /\t/)[0,1];
		$species = lc($species);
		$species =~ s/\s+$//;  $species =~ s/^\s+//;
		$display_name =~ s/\s+$//;  $display_name =~ s/^\s+//;
		push @lines,$species if($species);
		$name_disp{$species} = $display_name if($species && $display_name);
	}
	close(IN);

	my $test=$lines[0];
	###################################################################
	if ($test=~/^\d+/ && $test!~/[a-zA-Z\s]/) { # tax_id list	
		foreach  (@lines) {
			$id_name{$_}=0; #global
		}
		my $infile="$Bin/../dat/id_names";
		open (IN, $infile) || die ("can not find /dat/id_names\n");
		while (<IN>) {
			my $line=$_;
			if (/^(\d+)/) {
				my $id=$1;
				if (exists $id_name{$id}) { #global
					if ($line=~/\|\s([^\|\*]+)\s\*\sscientific name/) {
						$id_name{$id}=$1;   
					}
					
					
				}
			}
		}
		close(IN);

		print "\nTax_id\t\tTax_name\n\n" if (defined $Debug);
		#test whether there is unexsited id
		foreach  (sort {$a<=>$b} keys %id_name) {
			if (!$id_name{$_}) {
				$unfound{$_}=$id_name{$_};
				delete($id_name{$_});
			}else{
				print $_."\t\t".$id_name{$_}."\n" if(defined $Debug);
			}
		}
	}

	
	
	
	###################################################################
	if ($test=~/^[a-zA-Z]+/) {  # name list
		my %name_id; #local
		foreach  (@lines) {
			my $name=$_;
			$name=lc($name);
			$name_id{$name}=0; #local
		}
		my $infile="$Bin/../dat/id_names";
		open (IN,$infile ) || die ("can not find /dat/id_names\n");
		while (<IN>) {
			my $line=$_;
			$line=lc($line);
			my $id;
			if ($line=~/^(\d+)/) {
				$id=$1;
			}
			while ($line=~/\|\s+([^\|\*]+)\s+\*/g) {
				if ( exists $name_id{$1} ) {
					$name_id{$1}=$id;
					$id_name{$id}=$1; #global assignment
					last;
				}
			}
		}

		close(IN);
		
		print "\nTax_id\t\tTax_name\n\n" if (defined $Debug);
		#test whether there is unexsited id
		foreach  (keys %name_id) {
			if (!$name_id{$_}) {
				$unfound{$_}=$name_id{$_};
				delete($name_id{$_});
			}else{
				print $name_id{$_}."\t\t".$_."\n" if (defined $Debug);
			}
		}

	}


	
	###################################################################
	if ($test=~/^\d+\s+\w+/) {  # tax_id	  tax_name
		my %check;
		foreach  (@lines) {
			if (/(\d+)\s+(.+)/) {
				$id_name{$1}=$2; #globle
				$check{$1}=0;
			}
		}
		#check whether id exists
		my $infile="$Bin/../dat/id_names";
		open (IN, $infile) || die ("can not find /dat/id_names\n");
		while (<IN>) {
			my $line=$_;
			if (/^(\d+)/) {
				my $id=$1;
				if (exists $check{$id}) { #global
					if ($line=~/\|\s([^\|\*]+)\s\*\sscientific name/) {
						$check{$id}=$1;   
					}
					
					
				}
			}
		}
		close(IN);
		
		print "\nTax_id\t\tTax_name\n\n" if (defined $Debug);
		#test whether there is unexsited id
		foreach  (keys %check) {
			if (!$check{$_}) {
				$unfound{$_}=$id_name{$_};
				delete($id_name{$_});
			}else{
				print $_."\t\t".$id_name{$_}."\n" if(defined $Debug);
			}
		}
		
	}

}
####################################################


#construct %id_parents
####################################################
sub construct_id_parents{
	print "\nTax_id parent chain\n\n" if(defined $Debug);
	my $infile="$Bin/../dat/id_parents";
	open(IN,$infile) || die ("can't no do costruct_node\n");
	
	#first construct %id_parents #global
	while (<IN>) {
		my $line=$_;
		if (/(\d+)/) {
			if (exists $id_name{$1}) { #global
				my @temp=split(/\s+/,$line);
				@temp=reverse @temp;
				$id_parents{$1}=\@temp;
				if (defined $Debug) {
					print $1."\t|\t";
					foreach  (@temp) {
						print $_."\t";
					}
					print "\n";
				}
			}
		}
	}
	close(IN);
	
}
####################################################




#construct children links
####################################################
sub construct_node_link{

	foreach my $id (keys %id_parents) {
		my $pp=$id_parents{$id}; #global
		my $num=scalar(@$pp);
		for (my $i=0; $i<$num-1; $i++) {
			my $key="No".$id_parents{$id}[$i]."de"."link";
			my $val="No".$id_parents{$id}[$i+1]."de";
			if (!$node{$key}) { #global
				$node{$key}=[]; #node.link to a children list
			}
			
			my $pp=$node{$key}; #children list address
			my %temp; #store the children list
			foreach  (@$pp) {
				$temp{$_}=1;
			}

			if (!exists $temp{$val}) {
				push @$pp,$val;
			}
			
		}

	}
	


}
####################################################




#construct annotation links
####################################################
sub construct_node_anno{
	my %anno;   #which nodes need annotation
	foreach  (keys %id_parents) {
		my $pp=$id_parents{$_};
		foreach  (@$pp) {
			$anno{$_}=0;
		}
	}

	#search annotation of the specifide nodes which appear in the parent links
	my $infile="$Bin/../dat/id_information";
	open (IN,$infile) || die ("can not open ../dat/id_information\n");
	while (<IN>) {
		my $line=$_;
		if (/^(\d+)/) {
			my $id=$1;

			if (exists $anno{$id}) {
				
				my @ary = split(/\|/,$line);
				for (my $i=0; $i<@ary; $i++) {
					$ary[$i] =~ s/^\s+//;
					$ary[$i] =~ s/\s+$//;
				}
				my $var1=$ary[0];	$var1=~s/[^\w]/\_/g;
				my $var2=$ary[1];	$var2=~s/[^\w]/\_/g;
				my $var3=$ary[2];	$var3=~s/[^\w]/\_/g;
				$anno{$id}="[&&NHX".":T=".$var1.":S=".$var2.":R=".$var3."]";  #global assign
				
				$node{"No".$id."de"."anno"}=$anno{$id}; #close annotation
				

			}
		}
	}
	
	close(IN);
	

}
####################################################



#write the tree out in nhx format
####################################################
sub write_tree{

	print "\nThe process of tree constructing.....\n\n" if (defined $Debug);

	my $test=1; 
	while($test){
		my $len=length($tree); #$tree,global
		if(defined $Debug){
			for ( my $i=0; $i<$len; $i+=80 ) {
				my $str=substr($tree,$i,80);
				print $str."\n";
			}
			print "\n";
		}
		
		$test=0; #when all nodes are leaf nodes ,out of loop		
		my @inter;	
		while ($tree=~/(No\w+de)/g) { #pattern internal node 
			push @inter,$1;
		}

		foreach my $nodename (@inter) {
			if ($node{$nodename."link"} ne "") {
				my $pp=$node{$nodename."link"};
				my $content; #replace content
				
				foreach my $child (@$pp) {
					$content.=$child;
					$content.=$node{$child."anno"} if ($node{$child."anno"} ne "");
					$content.=",";

				}
				chop $content; #rm ","
				
				#judge whether this node is a leaf but have other leaf as children
				my $leaf=0; 
				if ($nodename=~/No(\d+)de/) {
					if (exists $id_name{$1}) { #globle
						$leaf=1;
						push @born,$1; #leaf nodes have chilren
					}
				}
				
				my $num=scalar(@$pp); #children number
				my $rep;
				if ($leaf==0 && $num==1) {
					$rep=$content;
					$tree=~s/$nodename\[[^\[\]]*\]/$rep/;
				}
				if ($leaf==0 && $num>1) {
					$rep="(".$content.")";
					$tree=~s/$nodename/$rep/;
				}
				if ($leaf==1) {
					$rep="(".$content.",".$empty.")";  #$empty is global
					$tree=~s/$nodename/$rep/;
				}
				
				$test=1;

			}
			
		}
	}



}
####################################################




####################################################
sub draw_tree{
	#replace tax_id with tax_name of each leaf node
	
	foreach  (keys %id_name) { #global
		my $tax_id=$_;
		my $tax_name=$id_name{$tax_id};
	
		$tax_name=~s/[^\w]/\_/g;
		$tree=~s/No($tax_id)de/$tax_name/;

	}

	$tree=~s/No0de/$id_name{0}/g;
	
	my $nhx = Tree::nhx_svg->new('is_real'=>0,show_inter=>1,width=>1000,skip=>20);
	$nhx->parse($tree);
	
	foreach my $p ($nhx->node) {
		$p->{S} =~ s/_/ /g;
		$p->{R} =~ s/_/ /g;
		if(defined $p->{C}){
			$p->{N} = $p->{S};
		}
		
		$p->{N} =~ s/_/ /g;
		$p->{N} = ucfirst($p->{N});
		
	}
	
	foreach my $p ($nhx->node) {
		$p->{N} = lc($p->{N});
		if (exists $name_disp{$p->{N}}) {
			$p->{N} = $name_disp{$p->{N}};
		}
		$p->{N} = ucfirst($p->{N});
	}
	
	mkdir($Outdir) unless(-d $Outdir);
	open OUT,">$Outdir/$file_core.nhx" || die "fail creat $Outdir/$file_core.nhx";
    print OUT $nhx->string($nhx->root,"nhx");
    close OUT;

	open OUT,">$Outdir/$file_core.nh" || die "fail creat $Outdir/$file_core.nh";
    print OUT $nhx->string($nhx->root,"nh");
    close OUT;

	open OUT,">$Outdir/$file_core.svg" || die "fail creat $Outdir/$file_core.svg";
    print OUT $nhx->plot;
    close OUT;
	

}
####################################################




#print out invalid inforamtion on screen
####################################################
sub disp_except{
	#print out unfound items 
	my $num=scalar(keys %unfound);  #global
	if ($num!=0) {
		print STDERR "\n\tUnfound id or name listed  and deleted:\n";
	}
	foreach  (keys %unfound) {
		print STDERR "\t\t".$_."\n";

	}

	#print out nodes which have child
	my $num=scalar(@born); #global
	if ($num!=0) {
		print STDERR "\n\tNodes has other node as children listed:\n";
	}
	foreach  (@born) {
		my $pp=$node{"No".$_."delink"};
		my $num=scalar(@$pp);
		print STDERR "\t".$_."\t".$id_name{$_}."\t"."have $num child\n";
		foreach  (@$pp) {
			my $id=$_;
			$id=~s/[Node]//g;
			print STDERR "\t\t".$id."\t".$id_name{$id}."\n";
		}
	}

}





__END__