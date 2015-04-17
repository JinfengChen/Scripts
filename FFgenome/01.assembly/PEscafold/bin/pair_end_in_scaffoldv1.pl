#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "./";
use SVG;
use FontSize;

my $font = FontSize->new();
my $font_family = "Arial";
my $font_size = 12;

if(@ARGV==0||($ARGV[0] ne "PRE"&&$ARGV[0] ne "DRAW"&&$ARGV[0] ne "ALL"))
{
   print "\tThis program is to draw the pair-ends relationship in the scaffold, you can process in one step like this:\n";
   print "\tperl $0 ALL <scaffold_list> <scaffold_seq> <soap_result_list>\n\n";
   print "\tor in two steps like this:\n";
   print "\tperl $0 PRE <scaffold_list> <scaffold_seq> <soap_result_list>\n";
   print "\tperl $0 DRAW\n\n";
   print "\tmodifying the last line in the file \"scaff_contig_pos.log\" can change the segment be drawn in the map\n";
   exit;
}

if($ARGV[0] eq "PRE"||$ARGV[0] eq "ALL")
{
   if(@ARGV<2||@ARGV>4)
   {
      print "perl $0 $ARGV[0] <scaffold_list> <scaffold_seq> <soap_result_list>\n";
      exit;
   }
   
   my %scaff_list=();
   open IN,"$ARGV[1]";
   while(<IN>)
   {
      chomp;
      $scaff_list{$_}=1;	
   }
   close IN;
   my $scaffold_seq=$ARGV[2];
   my $soap_re_list=$ARGV[3];
   &get_contig_pos($scaffold_seq,%scaff_list);
   &get_pair_info($soap_re_list,%scaff_list);
}

if($ARGV[0] eq "DRAW"||$ARGV[0] eq "ALL")
{
   &get_the_map();
}

sub get_the_map()
{
   my %hash_pair_info=();my %contig_compose=();my %segment_se=();my $name;my @tmp=();my $distance;
   open PAIR,"pair_ends_info.log", or die "Can't open the file pair_ends_info.log, maybe you should run the first";
   while(<PAIR>)
   {
      chomp;
      if($_=~/^\D/)
      {
      	$name=$_;
      }
      else
      {
         @tmp=split;
         $distance=$tmp[2]-$tmp[0];
         if($distance>=180&&$distance<420)
            {$hash_pair_info{$name}{300}{$tmp[0]}=$tmp[2];}
         elsif($distance>=420&&$distance<=980)
            {$hash_pair_info{$name}{700}{$tmp[0]}=$tmp[2];}
         elsif($distance>=1600&&$distance<=2400)
            {$hash_pair_info{$name}{2000}{$tmp[0]}=$tmp[2];}
         elsif($distance>=4000&&$distance<=6000)
            {$hash_pair_info{$name}{5000}{$tmp[0]}=$tmp[2];}
         elsif($distance>=8000&&$distance<=12000)
            {$hash_pair_info{$name}{10000}{$tmp[0]}=$tmp[2];}
         elsif($distance>=16000&&$distance<=24000)
            {$hash_pair_info{$name}{20000}{$tmp[0]}=$tmp[2];}
      }
   }
   close PAIR;
   
   my $flag=0;
   open CON,"scaff_contig_pos.log", or die "Can't open the file scaff_contig_pos.log, maybe you should run the first";
   while(<CON>)
   {
      chomp;
      if($_=~/^\D/)
      {
         if($flag==0)
            {$flag=1;}
         else
         {
            $segment_se{$name}{"start"}=$tmp[0];
            $segment_se{$name}{"end"}=$tmp[1];
         }
         $name=$_;
         $_=<CON>;
         @tmp=split(/\s+/,$_);
      }
      else
      {
         $contig_compose{$name}{$tmp[0]}=$tmp[1];
         @tmp=split;
      }
   }
   #$segment_se{$name}{"start"}=$tmp[0];
   #$segment_se{$name}{"end"}=$tmp[1];
   close CON;
   
   ################################################开始画图
   foreach my $name (keys(%contig_compose))
   {
      my $start=$segment_se{$name}{"start"};
      my $end=$segment_se{$name}{"end"};
      my $length=$end-$start;
      my $bp_per_pix=20;
      
      my $svg=SVG->new(width=>($length*3)/(9*$bp_per_pix),height=>500);
      my $path=$svg->group(
         style=>{fill=>'none','stroke-wide'=>2}
      );
      my $rect=$svg->group(
         style=>{stroke=>'black','stroke-wide'=>0}
      );
      
      ###########################################画出contig和坐标
      my $x_start=60;
      my $c_s=0;my $c_len=0;
      foreach my $key_1 (sort {$a<=>$b} keys %{$contig_compose{$name}})
      {
         if($key_1>=$start&&$key_1<=$end)
         {
         	$c_s=$key_1;
         	if($contig_compose{$name}{$key_1}<=$end)
         	   {$c_len=$contig_compose{$name}{$key_1}-$c_s;}
         	else
         	   {$c_len=$end-$c_s;}
         }
         elsif($key_1<$start&&$contig_compose{$name}{$key_1}>$start)
         {
         	$c_s=$start;
         	if($contig_compose{$name}{$key_1}<=$end)
         	   {$c_len=$contig_compose{$name}{$key_1}-$c_s;}
         	else
         	   {$c_len=$end-$c_s;}
         }
         
         $rect->rectangle(
            x=>$x_start+$c_s/$bp_per_pix,y=>450,
            width=>$c_len/$bp_per_pix,height=>10,
         ); 
      }
      &plot_ruler("svg",$svg,"Y",480, "X_start",$x_start,"X_end",$x_start+($end-$start)/$bp_per_pix,"bp_start",$start,"bp_end",$end,"scaletype","Kb","scaletypepos","right","scalestart","force","rulerstyle",2);
      
      ###########################################画出pair-end关系支持
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{20000}}))
      {
      	 my $len=($hash_pair_info{$name}{20000}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"black",}
         );
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{10000}}))
      {
      	 my $len=($hash_pair_info{$name}{10000}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"red",}
         );
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{5000}}))
      {
      	 my $len=($hash_pair_info{$name}{5000}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"blue",}
         );
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{2000}}))
      {
      	 my $len=($hash_pair_info{$name}{2000}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"green",}
         );
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{700}}))
      {
      	 my $len=($hash_pair_info{$name}{700}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"purple",}
         );
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$hash_pair_info{$name}{300}}))
      {
      	 my $len=($hash_pair_info{$name}{300}{$key_2}-$key_2)/$bp_per_pix;
      	 my $min=$x_start+$key_2/$bp_per_pix;
         $path->path(
            d=>"M$min 450 a1.2,1 0 1 1 $len,0",
            style=>{stroke=>"yellow",}
         );
      }
      open OUT,">$name\.svg";
      print OUT $svg->xmlify;
      close OUT;
      system("/home/jfchen/159/FFproject/tools/draw/svg2xxx_release/svg2xxx -t pdf $name.svg");
      system("/home/jfchen/159/FFproject/tools/draw/svg2xxx_release/svg2xxx -t png $name.svg");
   }
   
#   foreach my $key (keys(%segment_se))
#   {
#  	  print "$key\n";
#  	  foreach my $key2 (keys(%{$contig_compose{$key}}))
#  	  {
#  	  	 print "$key2\n";
#  	     print $contig_compose{$key}{$key2};
#  	  }
#  }

}

sub get_pair_info() #读取SOAP结果，得到pair-end关系的信息。
{
   my ($list,%pos)=@_;
   my %pair_info=();
   my $flag=0;
   my $allow_gap=50;
   
   open LIST,"$list";
   while(<LIST>)
   {
      chomp;
      open SOAP,"$_";
      while(<SOAP>)
      {
         chomp;
         my @tmp_1=split;
         if(defined($pos{$tmp_1[7]}))
         {
            my $line=<SOAP>;
            chomp($line);
            my @tmp_2=split(/\s+/,$line);

            my $distance=abs($tmp_1[8]-$tmp_2[8]);
            #########################################################################
            if($distance>=180&&$distance<420)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{300}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{300}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{300}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
            #########################################################################
            elsif($distance>=420&&$distance<=980)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{700}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{700}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{700}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
            #########################################################################
            elsif($distance>=1600&&$distance<=2400)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{2000}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{2000}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{2000}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
            #########################################################################
            elsif($distance>=4000&&$distance<=6000)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{5000}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{5000}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{5000}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
            #########################################################################
            elsif($distance>=8000&&$distance<=12000)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{10000}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{10000}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{10000}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
            #########################################################################
            elsif($distance>=16000&&$distance<=24000)
            {
               $flag=0;
               for(my $ii=-$allow_gap;$ii<=$allow_gap;$ii++)
               {
                  $flag=$flag||(defined($pair_info{$tmp_1[7]}{2000}{$tmp_2[8]+$ii}))
               }
               if(!$flag)
               {
                  if($tmp_1[8]>=$tmp_2[8])
                     {$pair_info{$tmp_1[7]}{2000}{$tmp_2[8]}=$tmp_1[8];}
                  else
                     {$pair_info{$tmp_1[7]}{2000}{$tmp_1[8]}=$tmp_2[8];}
               }
            }
               
            #####################################
            ############print "@tmp_1\n@tmp_2\n";
            #####################################
         }
      }
   }
   close LIST;
   
   open OUT,">pair_ends_info.log";
   foreach my $key_1 (keys(%pair_info))
   {
      print OUT "$key_1\n";
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{300}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{300}{$key_2}\n";
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{700}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{700}{$key_2}\n";
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{2000}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{2000}{$key_2}\n";
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{5000}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{5000}{$key_2}\n";
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{10000}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{10000}{$key_2}\n";
      }
      foreach my $key_2 (sort {$a<=>$b} keys (%{$pair_info{$key_1}{20000}}))
      {
         print OUT "$key_2 => $pair_info{$key_1}{20000}{$key_2}\n";
      }
   }
   close OUT;
}

sub get_contig_pos() 
{
   my $flag_1=0;my $flag_2=0;my $stat=0;my $pos=0;my $end=0;
   my ($fa,%pos)=@_;
   
   open OUT,">scaff_contig_pos.log";
   open IN,"$fa" or die "Can't open the file: $!\n";
   while(<IN>)
   {
   	  chomp;
   	  if($_=~/>(\w+)/)
   	  {
   	     if(defined($pos{$1})){
   	     	if($flag_1==1)
   	     	{
   	     	   if($flag_2==1){
                       print OUT "$pos\n";
                   }
   	     	   print OUT "1\t$pos\n";	
   	     	}
   	        $flag_1=1;$flag_2=0;$stat=0;$pos=0;
   	     }else{
                if($flag_1==1)
                {
                   if($flag_2==1){
                       print OUT "$pos\n";
                   }
                   print OUT "1\t$pos\n";
                }
                $flag_1=0;
             }
   	     print OUT "$1\n";
   	  }elsif($flag_1==1){
   	     my @seg=split(//,$_);
             for(my $i=0;$i<@seg;$i++){
                $pos++;
                if($seg[$i]=~/[ATCGatcg]/){
                   $flag_2=1;
                }else{
                   $flag_2=0;
                }
                if($flag_2!=$stat){
                  $end=$pos-1;
                  if($flag_2==1){
                     print OUT "$pos\t";
                  }else{
                     print OUT "$end\n";
                  }
                  $stat=$flag_2;
                }
             }
   	  }
     }
     #if($flag_2==1){
     #   print OUT "$pos\n";
     #   print OUT "1\t$pos\n";
     #}
     close IN;
}

sub plot_ruler()
{
        my %rulcfg = @_;
        $rulcfg{scaletype} ||= "scale";
        $rulcfg{scaletypepos} ||= "left";
        $rulcfg{scalestart} ||= "auto";
        $rulcfg{rulerstyle} ||= "3";


        my $scale_size = 6;
        my $divid = 50;
        my $unit_start;

        my $bp_len = $rulcfg{bp_end} - $rulcfg{bp_start};
        my $X_len = $rulcfg{X_end} - $rulcfg{X_start};

        my ($str,$str1,$str2,$unit);
        $str = $bp_len / $divid;
        $str = sprintf("%e",$str);
        if ($str=~/([0-9\.]+)e([0-9+-]+)/) {
                $str1 = $1;
                $str2 = $2;
        }
        $str1 = int ( $str1 + 0.5 );
        $unit = $str1 * 10 ** $str2;
        $unit = $rulcfg{bigscalesize}/10 if(defined $rulcfg{bigscalesize});

        my $g = $rulcfg{svg}->group('id'=>times().rand());

        ## draw the main axis
        $g->line('x1',$rulcfg{X_start},'y1',$rulcfg{Y},'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1);
        return if($rulcfg{bp_end} - $rulcfg{bp_start}  == 0);

        ##draw ruler mark text at the specified postion(left or right of the ruler)
        if ($rulcfg{scaletypepos} eq "left") {
                $g->text('x',$rulcfg{X_start}-$font->stringWidth($font_family,$font_size,$rulcfg{scaletype})-6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$font_family,"font-size",$font_size,"fill",'#000000');
        }
        if ($rulcfg{scaletypepos} eq "right") {
                $g->text('x',$rulcfg{X_end} + 6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$font_family,"font-size",$font_size,"fill",'#000000');
        }

        ##decide unit start
        if ($rulcfg{scalestart} eq "auto") {
                $unit_start = $rulcfg{bp_start} + ($unit - $rulcfg{bp_start} % $unit);
        }
        if ($rulcfg{scalestart} eq "force") {
                $unit_start = int($rulcfg{bp_start} / 10 + 0.5) * 10; ##....................?	
        }
        ## draw small scale lines
        for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit) {
                my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
                $g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '2');
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
        }
        ## draw big scale lines and text scales
        if((int($unit_start/1000))%10<5 && (int($unit_start/1000))%10>0)
        {
                                my $left=$unit_start%1000;
                                my $big=int$unit_start/10000;
                                $unit_start=$big*10000+5000+$left;
        }
        elsif((int($unit_start/1000))%10>5)
        {
                        my $left=$unit_start%1000;
                        my $big=int$unit_start/10000;
                        $unit_start=($big+1)*10000+$left;
        }
        else{$unit_start=$unit_start;}
        for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit*10) {
                my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
                $g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1)  if ($rulcfg{rulerstyle} eq '2');
                $g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
                my $scale_num = int ($i / 1000);
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$font->stringHeight($font_size),'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '1');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}+$scale_size+$font->stringHeight($font_size)+2,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '2');
                $rulcfg{svg}->text('x',$X - $font->stringWidth($font_family,$font_size,$scale_num) / 2,'y',$rulcfg{Y}-$scale_size-$font->stringHeight($font_size)+6,'fill','#000000','-cdata',$scale_num,'font-size',$font_size, 'font-family',$font_family) if ($rulcfg{rulerstyle} eq '3');
        }
}
