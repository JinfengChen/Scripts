#!/usr/bin/perl
## Change the file name in this script and run it in bin. 
## GFF file or BED file will be split into files for each chromosome in a chr directory.

###split rice gff file of gene anontation and TE annotation into chromosome files
my $genegff="/share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/RAP3.gff3.nr.gff";
my $tegff="/share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/IRGSP.build5.RepeatMasker.out.gff";
my $shootbed="/share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/IRGSPShoot.bed";
my $methy="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/IRGSP.methylation.bed";
my $h3k4="/home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OS.H3K4.bed";

#dealgff($genegff);
#dealgff($tegff);
#dealgff($shootbed);
#dealgff($methy);
dealgff($h3k4);

######################################
sub dealgff {
my ($gff)=@_;
chomp $gff;
my $dir=$gff.".chr";
`mkdir $dir` unless (-d $dir);
open IN, "$gff" or die "$!";
while(<IN>){
    next if ($_ eq "");
    next if ($_ =~/^##/);
    my @unit=split("\t", $_);
    $unit[0]=~/(\d+)/;
    $unit[0]=$1;
    my $chr;
    if (length $unit[0] == 2){
       $chr="chr".$unit[0];
    }else{
       $chr="chr0".$unit[0];
    }
    $unit[0]=$chr;
    my $line=join("\t",@unit);
    open OUT, ">>$dir/$chr"; 
         print OUT "$line"; 
    close OUT;
}
close IN;
}
