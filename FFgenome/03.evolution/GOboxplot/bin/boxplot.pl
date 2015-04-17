#!/usr/bin/perl
use Getopt::Long;
use Statistics::R;

GetOptions (\%opt,"class:s","table:s","project:s","help");


my $help=<<USAGE;
Prepare 4boxplot files for draw boxplot for a selected class of GO function, molecular_function/cellular_component/biological_process.
*.4boxplot: data matrix for drawing boxplot in R.
*.name: head of 4boxplot, need to check if the last column is TAB, delete this TAB.
*.boxplot.r: Run this scripts in R to generate a pdf figure file.

Run: perl boxplot.pl -class BP -table GO.result -project OB2OS_Ks
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my %GOclass=(
    "MF" => "molecular_function",
    "CC" => "cellular_component",
    "BP" => "biological_process"
    );
our %hash;
our %len;
our %median;
& parsetable($opt{table},$GOclass{$opt{class}});
#my $R=Statistics::R->new();
#$R->startR;

my @go=sort {$len{$b} <=> $len{$a}} keys %len;
my $tempnum=$hash{$go[0]};
#######################add NAs to these short array###
my @temparray=split(",",$tempnum);
my %gene4go;
foreach(keys %hash){
    my @array=split(",",$hash{$_});
    $gene4go{$_}=@array;  ### gene number for each GO term
    if (@array < @temparray){
       my $add=@temparray-@array;
       for(my $i=0;$i<$add;$i++){
          $hash{$_}.=",NA";
       }
    }
}

my $gonumber=0;
open OUT, ">$opt{class}.$opt{project}.GO.4boxplot" or die "$!";
#my $tmp=join("\t",@go);
#print OUT "$tmp\n";
#########sort go by median###############
@go=sort {$median{$b} <=> $median{$a}} keys %median;
foreach (@go){
  unless ($gene4go{$_} <= 40){
      print OUT "$_\t";
      $gonumber++;
  }
}
print OUT "\n";
for(my $i=0;$i<@temparray;$i++){
   for(my $j=0;$j<@go;$j++){
       if ($gene4go{$go[$j]} <= 40){
          next;
       }else{
          my @tmp=split(",",$hash{$go[$j]});
          print OUT "$tmp[$i]\t"
       }
   }
   print OUT "\n";
}
close OUT;
print "GO number\t$gonumber\n";
system ("head -n 1 $opt{class}.$opt{project}.GO.4boxplot > $opt{class}.$opt{project}.GO.4boxplot.name");


open OUT, ">$opt{class}.$opt{project}.GO.4boxplot.r" or die "$!";

print OUT "matrix(scan(\"$opt{class}.$opt{project}.GO.4boxplot\",skip=1),ncol=$gonumber,byrow=T) -> x","\n";
print OUT "scan(\"$opt{class}.$opt{project}.GO.4boxplot.name\",what=\"character\",sep=\"\\t\") -> name","\n";
print OUT 'colnames(x) <- name',"\n";
print OUT "pdf(\"$opt{class}.$opt{project}.GO.4boxplot.pdf\")","\n";
print OUT 'par(cex.axis=0.7,cex.main=0.7,cex.lab=0.7,omi=c(0,2,0,0),las=1)',"\n";
print OUT 'boxplot(as.data.frame(x),horizontal = TRUE,outline=FALSE,xlab="Ka/Ks (nonsynonymous/synonymous substitutions)")',"\n";
print OUT 'dev.off()',"\n";

close OUT;
############construct matrix#########################################
=pod
$R->send("x <- cbind(c($tempnum))");
$R->send("print (x)");
my $x=$R->read;
print "$x\n";
for(my $i=1;$i<@go;$i++){ 
    print "$go[$i]\n";
    $tempnum=$hash{$go[$i]};
    $R->send("x <- cbind(x,c($tempnum))"); 
}
$R->send("print (x)");
my $x=$R->read;
print "$x\n";
$R->send(q`pdf("/share/raid12/chenjinfeng/FFgenome/evolution/GOenrichTest/bin/GO.boxplot.pdf")`);
$R->send(q`boxplot(as.data.frame(x),main="GO.boxplot.horizontal",horizontal = TRUE)`);
$R->send(q`dev.off()`);
=cut

#####################sub functions#####################
sub parsetable{
my ($table,$class)=@_;
open IN, "$table" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t",$_);
   if ($unit[2] eq $class){
      my $id=$unit[3]."($unit[4])";
      $hash{$id}=$unit[7];
      $len{$id}=$unit[4];
      $median{$id}=$unit[6];
      #print "$id\t$hash{$id}\n";
      #print "$unit[3]\t$len{$unit[3]}\n";
   }
}
close IN;
}
