#!/usr/bin/perl

##  multi-process.pl -- control processes running with specified CPU number, 
##	read command from a file, each process can have more than one orderd
##	commands, but they should be put in the same one line, seperated by 
##	semicolon, for example: echo hello; sleep 3; perl test.pl; 

# Author: Fan Wei, fanw@genomics.org.cn

## Date: 2007-1-18

## 本程序的核心想法来源于李俊。
## 注意调用fork的时候非常危险，请万勿随意修改，本程序曾把北京的大型机搞死过两次。

use strict;
use Getopt::Long;

my $program_name=$1 if($0=~/([^\/]+)$/);
my $usage=<<USAGE; #******* Instruction of this program *********# 

Program: control multiple processes running

Usage: perl $program_name  <command_file>  
	-cpu <int>	number of CPU to use, default=3
	-cmd		output the commands but not execute
	-verbose	output information of running progress
	-help		output help information to screen

USAGE

my %opts;
GetOptions(\%opts, "cpu:s","cmd!","verbose!","help!");
die $usage if ( @ARGV==0 || defined($opts{"help"}));

#****************************************************************#
#--------------------Main-----Function-----Start-----------------#
#****************************************************************#

my $owner = `whoami`;

my @cmd;
while (<>) {
	chomp;
	s/&//g;
	next if(! /\S+/);
	push @cmd,$_;
}

if(exists $opts{cmd}){
	foreach  (@cmd) {
		print $_."\n";
	}
	exit;
}

Multiprocess(\@cmd,$opts{cpu});



#****************************************************************#
#------------------Children-----Functions-----Start--------------#
#****************************************************************#

##use mutiple CPUs at the same time
#########################################
sub Multiprocess{
	my $cmd_ap = shift;
	my $max_cpu = shift || 3; 
	my $total = @$cmd_ap;
	
	print STDERR "\n\tcmd num:  $total\n\tcpu num:  $max_cpu\n\n" if ( exists $opts{verbose} );

	my %report;
	for (my $i=1; $i<=10; $i++) {
		$report{int $total/10*$i}=$i*10;
	}
	
	for (my $i=0; $i<$total; $i++) {
		printf STDERR ("\tthrow out  %d\%\n",$report{$i+1}) if (exists $opts{verbose} && exists $report{$i+1});
		my $cmd=$$cmd_ap[$i];
		if ( fork() ) { 	
			wait if($i+1 >= $max_cpu); ## wait unitl all the child processes finished
		}else{          
			exec $cmd;  #child process
			exit();     #child process
		}
		sleep 0.1;
	}

	while (wait != -1) { sleep 1; }

	print STDERR "\tAll tasks done\n";
}
