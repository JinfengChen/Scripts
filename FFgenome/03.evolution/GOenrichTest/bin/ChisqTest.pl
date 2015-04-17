#!/usr/bin/perl
use Getopt::Long;
use Statistics::ChisqIndep;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Do Chi sequare Test for a file.
Example: 
name	#Test in FF	#Ref in FF	#Test in Rice	#Ref in Rice
apoplast	30	20	39	30
carbohydrate binding	60	68	143	138

Run: perl ChisqTest.pl fortest.txt
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
} 

my @obs;
open IN, "$ARGV[0]" or die "$!";;
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    next if ($unit[1]=~/\D+/);
    next if ($unit[1] eq "" or $unit[2] eq "" or $unit[3] eq "" or $unit[4] eq "");
    @obs=([$unit[1],$unit[2]],[$unit[3],$unit[4]]);
    #my @obs=([52,19],[39,3]);
    my $chi=new Statistics::ChisqIndep;
    $chi->load_data(\@obs);
    print "$_\t";
    #$chi->print_summary();  
    if ($chi->{valid}){
       print "$chi->{p_value}\n"; 
       #print "Rows: ", $chi->{rows}, "\n"; 
       #print "Columns: ", $chi->{cols}, "\n"; 
       #print "Degree of Freedom: ", $chi->{df}, "\n";
       #print "Total Count: ", $chi->{total}, "\n";
       #print "Chi-square Statistic: ", $chi->{chisq_statistic}, "\t";
       #print "p-value: ", $chi->{p_value}, "\n";
       #print "Warning: some of the cell counts might be too low.\n" if ($chi->{warning});
    }
}



