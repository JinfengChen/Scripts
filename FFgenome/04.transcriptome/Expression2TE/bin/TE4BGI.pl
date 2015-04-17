#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
Convert Dongying format to BGI gff format
perl $0 -gff TE.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $dbid;
open IN, "$opt{gff}" or die "$!";
open OUT, ">OBa.manual.TE.gff";
while(<IN>){
    chomp $_;
    print OUT "$_\n" and next if ($_=~/^#/);
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[5]=".";
    $unit[7]=".";
    my $target=$unit[8];
    $dbid++;
    if ($unit[2]=~/LINE/i){
       $unit[8]="ID=$dbid;$target;Class=LINE;";
    }elsif($unit[2]=~/SINE/i){
       $unit[8]="ID=$dbid;$target;Class=SINE;";
    }elsif($unit[2]=~/LTR/i){
       my $class="LTR";
       if ($unit[2]=~/copia/i){
          $class="LTR/Copia";
       }elsif($unit[2]=~/gypsy/i){
          $class="LTR/Gypsy";
       }
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/Stow/i){
       my $class="DNA/TcMar-Stowaway";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/Tour/i){
       my $class="DNA/Tourist";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/MITE/i){
       $unit[8]="ID=$dbid;$target;Class=DNA/MITE;";
    }elsif($unit[2]=~/Hel/i){
       my $class="RC/Helitron";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/Mariner/i or $unit[2]=~/Tc1/i){
       my $class="DNA/Mariner";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/Harbinger/i or $unit[2]=~/PIF/i){
       my $class="DNA/Harbinger";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/hAT/i){
       my $class="DNA/hAT";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/CACTA/i){
       my $class="DNA/En-Spm";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }elsif($unit[2]=~/Mu/i){
       my $class="DNA/MuDR";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }else{
       my $class="Unknown";
       $unit[8]="ID=$dbid;$target;Class=$class;";
    }
    $unit[2]="Transoposon";
    my $line=join("\t",@unit);
    print OUT "$line\n";
}
close IN;
close OUT;
