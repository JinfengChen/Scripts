#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
Convert Dongying format to Andrea gff format
perl $0 -gff TE.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $dbid;
open IN, "$opt{gff}" or die "$!";
open OUT, ">Obrachyantha.TE.gff";
while(<IN>){
    chomp $_;
    print OUT "$_\n" and next if ($_=~/^#/);
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[5]=".";
    $unit[7]=".";
    $dbid++;
    if ($unit[2]=~/LINE/i){
       $unit[2]="LINE_element";
       $unit[8]="ID=$dbid;rclass=I;subclass=I;rorder=LINE;superfamily=NA;family=NA";
    }elsif($unit[2]=~/SINE/i){
       $unit[2]="SINE_element";
       $unit[8]="ID=$dbid;rclass=I;subclass=I;rorder=SINE;superfamily=NA;family=NA";
    }elsif($unit[2]=~/LTR/i){
       $unit[2]="LTR_retrotransposon";
       my $super="NA";
       if ($unit[2]=~/copia/i){
          $super="Copia";
       }elsif($unit[2]=~/gypsy/i){
          $super="Gypsy";
       }
       $unit[8]="ID=$dbid;rclass=I;subclass=I;rorder=LTR;superfamily=$super;family=NA";
    }elsif($unit[2]=~/Stow/i){
       $unit[2]="MITE";
       $super="Stowaway";
       $unit[8]="ID=$dbid;rclass=II;subclass=III;rorder=MITE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/Tour/i){
       $unit[2]="MITE";
       $super="Tourist";
       $unit[8]="ID=$dbid;rclass=II;subclass=III;rorder=MITE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/MITE/i){
       $unit[2]="MITE";
       $unit[8]="ID=$dbid;rclass=II;subclass=III;rorder=MITE;superfamily=NA;family=NA";
    }elsif($unit[2]=~/Hel/i){
       $unit[2]="helitron";
       $unit[8]="ID=$dbid;rclass=II;subclass=II;rorder=Helitron;superfamily=Helitron;family=NA";
    }elsif($unit[2]=~/Mariner/i or $unit[2]=~/Tc1/i){
       $unit[2]="terminal_inverted_repeat_element";
       $super="Mariner";
       $unit[8]="ID=$dbid;rclass=II;subclass=I;rorder=DNA_TE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/Harbinger/i or $unit[2]=~/PIF/i){
       $unit[2]="terminal_inverted_repeat_element";
       $super="Harbinger";
       $unit[8]="ID=$dbid;rclass=II;subclass=I;rorder=DNA_TE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/hAT/i){
       $unit[2]="terminal_inverted_repeat_element";
       $super="hAT";
       $unit[8]="ID=$dbid;rclass=II;subclass=I;rorder=DNA_TE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/CACTA/i){
       $unit[2]="terminal_inverted_repeat_element";
       $super="CACTA";
       $unit[8]="ID=$dbid;rclass=II;subclass=I;rorder=DNA_TE;superfamily=$super;family=NA";
    }elsif($unit[2]=~/Mu/i){
       $unit[2]="terminal_inverted_repeat_element";
       $super="Mutator";
       $unit[8]="ID=$dbid;rclass=II;subclass=I;rorder=DNA_TE;superfamily=$super;family=NA";
    }else{
       $unit[2]="dispersed_repeat";
       $unit[8]="ID=$dbid;rclass=UNKNOWN;subclass=UNKNOWN;rorder=UNKNOWN;superfamily=UNKNOWN;family=UNKNOWN";  
    }
    my $line=join("\t",@unit);
    print OUT "$line\n";
}
close IN;
close OUT;
