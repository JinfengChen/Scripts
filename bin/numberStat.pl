#!/usr/bin/perl

## numberStat.pl -- do basic statistics on a set of numbers

## this program is designed for common number statistics, such as:
## count number and caculate total vaule, caculate mean and median value,
## caculate Standard Deviation(SD), find the minimum and maximum vaule,
## caculte N50 vaule (this is only used for sequence assembly area)

## Note that the program is designed to work through the freely combination
## of options, eg, "-type mean:max" will caculate the mean value and maximum
## value; by default you do not need to set the type option, then it will 
## caculate all types of values, and give a detailed result report.

# Author: Fan Wei, fanw@genomics.org.cn

## Date: 2007-1-18

#Include Modules
use strict;
use Getopt::Long; 

#Instructions and advisions of using this program
my $program=`basename $0`;
chomp $program;
my $usage=<<USAGE;

Program: do basic statistics on a set of numbers

Usage: $program  <STDIN or infile>  
	-type		mean:median:SD:max:min:number:total:N50:all, default=all
	-help		output help information to screen

Examples: 
	cat num_file.txt | $program; 
	cat num_file.txt | $program -type median
	cat num_file.txt | $program -type mean:median:max

USAGE

#################-Main-Function-#################

#get options and parameters
my %opts;
GetOptions(\%opts,"type:s","help!");
die $usage if ( exists $opts{"help"} );
$opts{type} = exists $opts{type} ? $opts{type} : "all";

##Constant and global variables
my (@num,$number,$total,$mean,$median,$SD,$add_square,$max,$min,$N50);

## read numbers from infile or pipe
while (<>) {
	chomp;  ## 去掉换行符
	next if(! /\S+/); ## 忽略空白行
	$_ = $1 if(/^(\S+)/); ## 取每行的第一个字段
	next if(/[^-\d\.eE]/); ## 忽略非正常数字的错误行
	push @num,$_;
}

@num=sort {$a<=>$b} @num;
my $loop=0;
foreach  (@num) {
	$total+=$_;
	$add_square+=$_*$_;
	$loop++;
}
$number=@num;
$min=$num[0];
$max=$num[$number-1];
$mean=$total/$number;
$median=$num[int $number/2];
$SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );

## caculate N50 value
my $add_value = 0;
my $half_total = $total / 2;
for (my $i=@num-1; $i>=0; $i--) {
	$add_value += $num[$i];
	if ($add_value >= $half_total) {
		$N50 = $num[$i];
		last;
	}
}

##output the result
if ($opts{type} eq "all") {
	print "total number:  $number\n";
	print "total value:   $total\n";
	print "mean value:    $mean\n";
	print "median value:  $median\n";
	print "SD value:      $SD\n";
	print "max value:     $max\n";
	print "min value:     $min\n";
	print "N50 value:     $N50\n";
	exit;

}else{
	my @key=split(/:/,$opts{type});
	my %hash=(mean=>$mean,median=>$median,SD=>$SD,number=>$number,total=>$total,max=>$max,min=>$min,N50=>$N50);
	foreach  (@key) {
		print $hash{$_}."\t" if(exists $hash{$_});
	}
	print "\n";
	exit;
}



#################-Sub--Routines-#################
