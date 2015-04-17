#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blasttable:s","help");


my $help=<<USAGE;
perl $0 --blasttable

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $ref=blastpair($opt{blasttable});
#my $ref=getblast($opt{blasttable});
print "QueryID\tHitID\tQueryCoverLength\tQueryCoverRate\tIdentity\tQueryAnnotation\n";
foreach (sort keys %$ref){
    print "$_\t$ref->{$_}->[0]\t$ref->{$_}->[1]\t$ref->{$_}->[2]\t$ref->{$_}->[3]\t$ref->{$_}->[4]\n";
    #print "$_\t$ref->{$_}\n";
}

sub getblast
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    next if ($unit[13] > 1e-10);
    if ($unit[0]=~/\w+/){
       my $len=$unit[1];
       my $hitlen=$unit[3]-$unit[2]+1;
       my $coverage=$hitlen/$len;
       my $hit=$unit[4];
       my $identity=$unit[8];
       my $subjectlen=$unit[5];
       $hash{$unit[0]}="$hit\t$hitlen\t$coverage\t$identity" unless exists $hash{$unit[0]};
    }
}
close IN;
return \%hash;
} 


sub blastpair
{
####read a blasttable file, return a ref of hash contain pair of homologous gene that have e-value <=1e-15,
####and coverage of > 70% for both protein sizes;
my ($blast)=@_;
my %record;
open IN, "$blast" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   if ($unit[13] > 1e-5){
       next;
   }
   next if ($unit[0] eq $unit[4]); ## remove self match
   my $pair=$unit[0].":".$unit[4];
   if (exists $record{$pair}){
        my $refarray=$record{$pair};
        push (@$refarray,[@unit]);
        $record{$pair}=$refarray;
   }else{
        my @array;
        push (@array,[@unit]);
        $record{$pair}=\@array;
   }
}
close IN;

my %homologous;
foreach (keys %record){
      #print "$_\n";
      #print Dumper($record{$_}),"\n";
      my $refarray=$record{$_};
      my $num=@$refarray;
      my $hsp=1;
      my $queryid =$$refarray[0][0];
      my $hitid   =$$refarray[0][4];
      my $bitscore=$$refarray[0][12];
      my $queryanno=$$refarray[0][14];
      my $identity=$$refarray[0][8];
      my $querylen=$$refarray[0][1];
      my $hitlen  =$$refarray[0][5];
      my $longestq=$$refarray[0][3]-$$refarray[0][2]+1;
      my $longesth=$$refarray[0][7]-$$refarray[0][6]+1;
      my $matchq  =$$refarray[0][3]-$$refarray[0][2]+1;
      my $matchh  =$$refarray[0][7]-$$refarray[0][6]+1;
      if ($num > 1){
          my $hspstartq=$$refarray[0][2];
          my $hspstarth=$$refarray[0][6];
          my $hspendq  =$$refarray[0][3];
          my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
                      next;
                }else{
                      $hsp++;
                      $hspendq=$$refarray[$i][3];
                      $hspendh=$$refarray[$i][7];
                      $bitscore+=$$refarray[$i][12];
                      $identity+=$$refarray[$i][8];
                      $matchq  +=$$refarray[$i][3]-$$refarray[$i][2]+1;
                      $matchh  +=$$refarray[$i][7]-$$refarray[$i][6]+1;
                }
          }
          $longestq=$hspendq-$hspstartq+1;
          $longesth=$hspendh-$hspstarth+1;
          $identity=$identity/$hsp;
      }
      my $qhspcoverage=$matchq/$querylen;
      my $hhspcoverage=$matchh/$hitlen;
      my $qlencoverage=$longestq/$querylen;
      my $hlencoverage=$longesth/$hitlen;
      #if ($qhspcoverage >= 0.7 and $hhspcoverage >= 0.7){
      #if ($qlencoverage >= 0.7 and $hlencoverage >= 0.7){
      #if ($qhspcoverage >= 0.3 and $hhspcoverage >= 0.3 and $qlencoverage >= 0.5 and $hlencoverage >= 0.5){
      #    my $temp="$queryid\t$hitid";
          #print "$temp\n";
      #    $homologous{$temp}=$identity;
      #}
      ## if the match hsp length > this value in hash
      #if ($matchh > $homologous{$queryid}->[1] or $longesth > $homologous{$queryid}->[3]){
         #$homologous{$queryid}=[$hitid,$matchh,$hhspcoverage,$longesth,$hlencoverage];
      #if ($matchh > $homologous{$queryid}->[1] or $hhspcoverage > $homologous{$queryid}->[2]){
      #next if ($identity < 0.9);
      if ($matchh > $homologous{$queryid}->[1]){
         $homologous{$queryid}=[$hitid,$matchq,$qhspcoverage,$identity,$queryanno];
      }
}
return \%homologous;
}
 
