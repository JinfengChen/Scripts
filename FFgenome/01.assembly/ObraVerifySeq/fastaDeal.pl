#!/usr/bin/perl

=head1 Program Description

This program is designed to be a universal fasta format dealer, including attribution
statistics, large file cutting, word patterning, fragment substructing, and sequence
reforming.

For fasta file cutting, you can either specify the sequence number in each 
sub file or sepecify the sub file number you wanted; The number of files in each direcotry
can't exceed the "file_in_dir" value, if it exceeds, then the result files 
will be put into 2-rank directory. You can set the "file_in_dir" value as you want, default
is no limited.

For this group of options( --substruct --reverse  --complement ), it can only deal with
single sequence fasta format. If the input is multi-seqeunce fasta file, it will only recognize
the first seqeunce.


=head1 Contact & Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 4.1,  Date: 2007-8-5

=head1 Command-line Option
 
  --attribute <str>       get differnt attributions, id:head:seq:len:lenwogap:gc:x:lc:uc.   
  
  --cuts <num>            cut file with the specified number of seqeunces in each sub file 
  --cutf <num>            cut file with the specified number of subfile in total
  --file_in_dir <num>     the max file number in each directory, default is no limited
  --outdir <str>          set the result dir, only used for cuts and cutf, default=./

  --get_id <str>          get one sequence by specified ID
  --pattern <str>         get out sequence which can match the specifed word
  --unpattern <str>       get out sequence which do not match the specifed word

  --substruct <str>       100-200(for exa),get a fragment from the sequence based on start/end coordinates
  --reverse               reverse the sequence 
  --complement            complement the sequence

  --sample <num-num>      get sequences by specified order range. Only give one num will only extract one seqeunce. 
  --reform <str>          lowerize|upperize|line50|pure, reform/modify the sequence as you want
  
  --verbose               output verbose information to screen  
  --help                  output help information to screen  

=head1 Usage Exmples

  fastaDeal.pl -attr id:len:gc chrY.fa

  fastaDeal.pl --cuts 2 cds.fa 
  fastaDeal.pl --cutf 20 cds.fa

  fastaDeal.pl --pat "AK" cds.fa 
  fastaDeal.pl --unpat "AK" cds.fa

  fastaDeal.pl --sample 5-10 cds.fa
  fastaDeal.pl --sample 10 cds.fa

  fastaDeal.pl -sub 100-200 chr1.fa
  fastaDeal.pl -sub 100-200 -rev -com chr1.fa

  fastaDeal.pl -reform upperize  cds.fa 

=cut


use strict;
use Getopt::Long;
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my ($Attribute,$Cuts,$Cutf,$File_in_dir_num,$Pattern,$Unpattern,$Sample,$Get_id);
my ($Substruct,$Reverse,$Complement,$Reform,$Outdir);
my ($Verbose,$Help);
GetOptions(
	"attribute:s"=>\$Attribute,
	"cuts:n"=>\$Cuts,
	"cutf:n"=>\$Cutf,
	"outdir:s"=>\$Outdir,
	"file_in_dir:n"=>\$File_in_dir_num,
	"get_id:s"=>\$Get_id,
	"pattern:s"=>\$Pattern,
	"unpattern:s"=>\$Unpattern,
	"sample:s"=>\$Sample,
	"substruct:s"=>\$Substruct,
	"reverse"=>\$Reverse,
	"complement"=>\$Complement,
	"reform:s"=>\$Reform,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
##$File_in_dir_num ||= 1000;

die `pod2text $0` if ($Help);

get_attribute() if($Attribute);

cut_fasta($ARGV[0],$File_in_dir_num) if($Cuts);

cut_fasta_advanced($ARGV[0]) if($Cutf); ##only one rank directory

pattern_head() if($Pattern || $Unpattern || $Get_id);

sample_seq() if($Sample);

substruct_frag() if($Substruct || $Reverse || $Complement);

reform_seq() if($Reform);

####################################################
################### Sub Routines ###################
####################################################

####################################################
sub reform_seq {
	
	$/=">";<>;$/="\n";
	while (<>) {
		my $title = $_;
		my $seq_name = $1 if($title =~ /^(\S+)/);
		
		$/=">";
		my $seq=<>;
		chomp $seq;
		$/="\n";
		
		print STDERR "length of $seq_name is 0\n" if($Verbose);

		$seq = lc($seq) if($Reform =~ "lowerize");
		$seq = uc($seq) if($Reform =~ "upperize");
		
		Display_seq(\$seq,$1) if($Reform =~ /line(\d+)/);
		
		if($Reform eq "pure"){
			$seq =~ s/[^a-zA-Z]//g;
			Display_seq(\$seq);
		}

		print ">".$title.$seq;
	}
}

sub substruct_frag{
	my ($seq_name,$head,$seq);
	
	$/=">"; <>; $/="\n";
	while (<>) {
		chomp;
		$head = $_;
		$seq_name = $1 if($head =~ /^(\S+)/);
		
		$/=">";
		$seq = <>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		last;
	}
	
	my $seq_len = length($seq);
	my $sub_head = ">".$seq_name;
	my $sub_seq = $seq;
	if($Substruct =~ /(start|\d+)-(end|\d+)/){
		my ($sub_start,$sub_end) = ($1,$2);
		$sub_start = 1 if($sub_start eq 'start' || $sub_start < 1);
		$sub_end = $seq_len if($sub_end eq 'end' || $sub_end > $seq_len);	
		$sub_head .= "_".$sub_start."_".$sub_end;
		$sub_seq = substr($seq,$sub_start-1,$sub_end-$sub_start+1);
	}
	
	$sub_seq = reverse $sub_seq if($Reverse);
	$sub_head .= "_reverse" if($Reverse);
	$sub_seq =~ tr/AGCTagct/TCGAtcga/ if($Complement);
	$sub_head .= "_complement" if($Complement);
	
	Display_seq(\$sub_seq);
	
	print $sub_head."\n".$sub_seq;
}

#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################



####################################################
sub sample_seq{
	
	my ($start,$end);
	if ($Sample =~ /(\d+)-(\d+)/){
		($start,$end) = ($1,$2); 
	}elsif ($Sample =~ /(\d+)/){
		($start,$end) = ($1,$1); 
	}else{
		die "please use --help to get help\n";
	}

	my $loop;
	$/=">";<>;$/="\n";
	while (<>) {
		$loop++;
		my $title = $_;
		
		$/=">";
		my $seq=<>;
		chomp $seq;
		$/="\n";

		print ">".$title.$seq if($loop >= $start);
		last if($loop >= $end);
	}
	
}


####################################################
sub pattern_head{
	
	$/=">";<>;$/="\n";
	while (<>) {
		my $title=$_;  
		chomp $title;
		
		$/=">";
		my $seq=<>;
		chomp $seq;
		$/="\n";
		chomp $seq;
		
		if ($Pattern) {
			my $pat=$Pattern;
			if ($pat=~/^\^/) {
				$pat=~s/^\^//;
				print ">".$title."\n".$seq."\n"  if($title=~/^$pat/);
			}elsif($pat=~/\$$/){
				$pat=~s/\$$//;
				print ">".$title."\n".$seq."\n"  if($title=~/$pat$/);
			}else{
				print ">".$title."\n".$seq."\n"  if($title=~/$pat/);
			}

		}
		if ($Unpattern) {
			my $pat=$Unpattern;
			if ($pat=~/^\^/) {
				$pat=~s/^\^//;
				print ">".$title."\n".$seq."\n"  if($title!~/^$pat/);
			}elsif($pat=~/\$$/){
				$pat=~s/\$$//;
				print ">".$title."\n".$seq."\n"  if($title!~/$pat$/);
			}else{
				print ">".$title."\n".$seq."\n"  if($title!~/$pat/);
			}
		}

		if ($Get_id) {
			my $seq_id = $1 if($title =~ /^(\S+)/);
			if ($Get_id eq $seq_id) {
				print ">".$title."\n".$seq."\n";
				last;
			}
		}
	}

}



##get different attributions of each sequence 
##id:head:seq:len:lenwogap:gc:x:lc:uc
##lenwogap, lenght exclude gap, suppose N n as gap
##gc, includes g c G c
##x, includes x X
##lc, includes a g c t
##uc, includes A G C T
####################################################
sub get_attribute{
	
	my @attr = split(/:/,$Attribute);
	my (@key,%work);
	foreach  (@attr) {
		my $str = uc($_); ## upperize all the attribution words
		$work{$str} = 1;
		push @key,$str;
	}
	
	$/=">";<>;$/="\n";
	while (<>) {
		chomp;
		my %attr;		
		$attr{HEAD} = $_ if(exists $work{HEAD});
		$attr{ID} = $1 if(exists $work{ID} && /^(\S+)/);
		
		$/=">";
		my $seq=<>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		my $len=length($seq);
		my $N = $seq=~tr/Nn//; ## stand for gap or low quality
		my $Lenwogap = $len - $N;

		$attr{SEQ}=$seq if(exists $work{SEQ});
		$attr{LEN}=$len if(exists $work{LEN});
		$attr{LENWOGAP}= $Lenwogap if(exists $work{LENWOGAP});


		if (exists $work{GC}) {
			my $GC = $seq=~tr/GCgc//; ## stand for GC content
			$attr{GC} = ($Lenwogap>0) ? ($GC / $Lenwogap) : "none";
		}
		
		if (exists $work{X}) {
			my $X = $seq=~tr/Xx//; ## stand for masked region
			$attr{X} = ($Lenwogap>0) ? ($X / $Lenwogap) : "none";
		}
		
		if (exists $work{LC}) {
			my $lc = $seq=~tr/agct//; ## stand for masked region
			$attr{LC} = ($Lenwogap>0) ? ($lc / $Lenwogap) : "none";
		}
		if (exists $work{UC}) {
			my $uc = $seq=~tr/AGCT//; ## stand for masked region
			$attr{UC} = ($Lenwogap>0) ? ($uc / $Lenwogap) : "none";
		}


		## make result in specified order
		my $output; 
		foreach  (@key) {
			$output .= $attr{$_}."\t";
		}
		chop $output;
		$output .= "\n";
		print $output;
	}

}

#cut fasta file, can't accept STDIN
####################################################
sub cut_fasta_advanced{
	my $infile = shift;
	
	##calculate the total sequence length in a fasta file
	my $Total_len=0;
	open IN,$infile or die "can not open $infile:$!";
	$/=">";<IN>;$/="\n";
	while(<IN>) {
			chomp;
			my $head=$_;
			$/=">";
			my $seq=<IN>;
			chomp $seq;
			$/="\n";
			$seq=~s/\s//g;
			$Total_len+=length($seq);

	}
	close IN;
	##warn $Total_len;
	my $file_name = `basename $infile`; ## only file name of $infile
	chomp $file_name;
	my $out_dir = (defined $Outdir) ? "$Outdir/$file_name.cut" : "./$file_name.cut";
	`rm -r $out_dir` if(-e "$out_dir");
	mkpath("$out_dir");

	my $Sub_len=int ($Total_len/$Cutf);
	##warn $Sub_len;
	my $Cur_len=0;
	my $Cur_content;
	my $file_mark = 1;
	open IN,$infile or die "can not open $infile:$!";
	$/=">";<IN>;$/="\n";
	while(<IN>) {
			chomp;
			my $head=$_;
			$/=">";
			my $seq=<IN>;
			chomp $seq;
			$/="\n";
			
			my $str = $seq;
			$str=~s/\s//g;
			$Cur_len+=length($str);
			$Cur_content .= ">$head\n$seq";

			if ($Cur_len >= $Sub_len) {
				open OUT,">$out_dir/$file_name.$file_mark" || die "fail";
				print OUT $Cur_content;
				close OUT;
				$Cur_content = "";
				$Cur_len = 0;
				$file_mark++;
			}
	}
	
	##make the last file
	if ($Cur_content) {
		open OUT,">$out_dir/$file_name.$file_mark" || die "fail";
		print OUT $Cur_content;
		close OUT;
	}

}

#cut fasta file, can't accept STDIN
####################################################
sub cut_fasta{
	my $infile = shift;
	my $file_in_dir_num = shift;
	
	my $file_name = `basename $infile`; ## only file name of $infile
	chomp $file_name;
	my $out_dir = (defined $Outdir) ? "$Outdir/$file_name.cut" : "./$file_name.cut";
	my $total_num = `grep -c '>' $infile`;  chomp $total_num;

	my ($seq_num,$file_num);
	if ($Cuts) {
		$seq_num = $Cuts;
		$file_num = int($total_num / $seq_num) + 1;
	}
	
	my $dir_rank = (defined $file_in_dir_num && $file_num > $file_in_dir_num) ? 2 : 1;
	
	print STDERR  "\nseq number in each sub file:   $seq_num\n".
			"sub file number in total:      $file_num\n".
			"sub directory rank:            $dir_rank\n\n\n" if($Verbose);

	`rm -r $out_dir` if(-e "$out_dir");
	mkpath("$out_dir");
	
	#my ($file_name,$suffix) = ($1,$2) if ($infile=~/(\S+)\.([^\.]+)$/);
	

	my $mark=$file_num; ## sub file numeric mark
	$mark=~tr/[0-9]/0/;
	$mark++;
	my $submark;
	$submark = int($file_num / $file_in_dir_num) + 1 if(defined $file_in_dir_num);
	$submark=~tr/[0-9]/0/; ## sub dir numeric mark
	
	my $num_loop=0;
	open IN, $infile || die "$infile can be open\n";
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $title=">".$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";

		if ($dir_rank == 2 && $num_loop % ($seq_num * $file_in_dir_num) == 0) {
			$submark++;
			mkdir("$out_dir/$submark");
		}
		
		if($num_loop % $seq_num == 0){
			print STDERR "generating:  $out_dir/$file_name.$mark\n" if( $Verbose && $dir_rank == 1);
			open ( OUT,">$out_dir/$file_name.$mark" ) if($dir_rank == 1);

			print STDERR "generating:  $out_dir/$submark/$file_name.$mark\n" if($Verbose && $dir_rank == 2);
			open ( OUT,">$out_dir/$submark/$file_name.$mark" ) if($dir_rank == 2);
			
			$mark++;
		}

		print OUT $title.$seq;
		close ( OUT ) if( ($num_loop % $seq_num) == ($seq_num-1) );
		
		$num_loop++;
	}

}
