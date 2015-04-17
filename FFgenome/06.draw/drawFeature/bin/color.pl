#!/usr/bin/perl

use SVG;
my %ref=(
       "Gene" => [0,20],
       "DNA"  => [0,25],
       "LTR"  => [0,40],
       "Other"=> [0,10]
);
my %qry=(
       "Gene" => [0,20],
       "DNA"  => [0,25],
       "LTR"  => [0,40],
       "Other"=> [0,10]
);
my $refhash=ref2array(\%ref);
my $qryhash=ref2array(\%qry);
my @color=colortab();

my $svg=SVG->new(width=>800,height=>400);
my $y=50;
my $x=30;
for(my $i=0;$i<26;$i++){
       my $n=int($i*1018/26);
       #print "$color[$n]\n";
       my ($red,$yellow,$blue)=split("\t",$color[$n]);
       print "$n\t$red\t$yellow\t$blue\n";
       my $tag=$svg->rectangle(
                    x=>$x,y=>$y,
                    width=>20,height=>10,
                    style=>{
                        fill=> "rgb($red,$yellow,$blue)"
                    }
       );
       my $refx1=$x;
       foreach (keys %$refhash){ 
         $refx1+=80;
         my $refy1=$y+8;
         my $genehash=$refhash->{$_};
         my $geneint =$genehash->{$i};
         my $head=$svg->text(
                    x=>$refx1,y=>40
         )->cdata($_);
         my $text=$svg->text(
                    x=>$refx1,y=>$refy1
         )->cdata($geneint);
       }
       my $note1=$svg->text(
                x=>220,y=>20
       )->cdata(Japonica);
       my $qryx1=$refx1+20;
       foreach (keys %$qryhash){
         $qryx1+=80;
         my $qryy1=$y+8;
         my $genehash=$qryhash->{$_};
         my $geneint =$genehash->{$i};
         my $head=$svg->text(
                    x=>$qryx1,y=>40
         )->cdata($_);
         my $text=$svg->text(
                    x=>$qryx1,y=>$qryy1
         )->cdata($geneint);
       }
       my $note2=$svg->text(
                x=>550,y=>20
       )->cdata(Brachyantha);
       if ($y + 10 >= 600){
          $y=100;
          $x=$x+100;
       }else{
          $y+=10;
       }
}
my $annotation="Proportion of basepair in 1 Mb windows with 200kb steps along chromosome";
my $note3=$svg->text(
     x=>110,y=>330
)->cdata($annotation);

writesvg("test.svg",$svg);


sub ref2array {
my ($ref)=@_;
my %refhash;
foreach (keys %$ref){
  my %hash;
  my $min=$ref->{$_}->[0];
  my $max=$ref->{$_}->[1];
  print "$_\t$min\t$max\n";
  my $counter=0;
  my $interval=($max-$min)/26;
  for (my $i=$min;$i<$max;$i+=$interval){
        my $limit=sprintf ("%.2f",$i+$interval);
        $i2=sprintf ("%.2f",$i);
        print "$counter\t$i2\t$limit\n";
        $hash{$counter}="$i2"."-"."$limit";
        $counter++;
  }
  $refhash{$_}=\%hash; 
}
return \%refhash;
}
#writesvg("test.svg",$svg);

sub colortab
{
        my ( @color, @arr );
        @arr= (0,0,128);
        for ( my $i=128; $i<255; $i++ ){
                push (@color, "$arr[0]\t$arr[1]\t$arr[2]");
                $arr[2]++;
        }
        for ( my $i=0; $i<255; $i++ ){
                $arr[1]++;
                push (@color, "$arr[0]\t$arr[1]\t$arr[2]");
        }
        for ( my $i=0; $i<255; $i++ ){
                $arr[0]++;
                $arr[2]--;
                push (@color, "$arr[0]\t$arr[1]\t$arr[2]");
        }
        for ( my $i=0; $i<255; $i++ ){
                $arr[1]--;
                push (@color, "$arr[0]\t$arr[1]\t$arr[2]");
        }
        for ( my $i=0; $i<127; $i++ ){
                $arr[0]--;
                push (@color, "$arr[0]\t$arr[1]\t$arr[2]");
        }

return @color;
}


sub writesvg {
my ($file,$svg)=@_;
#print "$file\n";
open OUT, ">$file" or die "can not open my file";
       print OUT $svg->xmlify;
close OUT;
       system "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx $file -m 400 -t pdf";
}
