#!/usr/bin/perl
## Based on the relation from joinScaffold.txt, we check the axt file and move out these scaffold that do not belong to chr

my $join="/share/raid12/chenjinfeng/FFgenome/assmbly/superscaf/input/joinScaffold.txt";


my %scaf2chr;
open IN, "$join" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_); 
   my $scaf=sprintf("Scaffold%06d",$unit[1]);
   $unit[0]=~s/\D+//g;
   my $chr =sprintf("chr%02d",$unit[0]);
   unless (exists $scaf2chr{$scaf}){
      $scaf2chr{$scaf}=$chr;
   }
}
close IN;

my @axt=glob("./*.axt");
foreach(@axt){
   my $file=$_;
   my $out=$file.".check";
   if ($_=~/(\w+).axt.chain.prenet.net.filter.net.axt/){
      my $chr=$1;
      open FILE, "$file" or die "$!";
      open OUT, ">$out" or die "$!";
             while (<FILE>){
                my @word;
                if ($_=~/^\d+/){
                  @word=split(" ",$_);      
                  if ($scaf2chr{$word[4]} ne $chr or $word[7] ne "+"){
                        <FILE>;
                        <FILE>;
                        <FILE>;
                  }else{
                        print OUT "$_";
                        for(my $i=0;$i<=2;$i++){
                           my $line=<FILE>;
                           print OUT "$line";
                        }
                  }
                }
             }
      close OUT;
      close FILE; 

   } 

}


