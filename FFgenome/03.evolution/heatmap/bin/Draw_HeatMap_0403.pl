#!/usr/bin/perl
# usage: perl $0 A01.data.distri A01.data.densi &
#use strict;
use Data::Dumper;
use lib "$Bin/../lib";
use SVG;

my ( $newwidth, $newheight, $edge_dis ) = ( 1200, 400, 50 );
my $svg=SVG->new( 'width',$newwidth,'height',$newheight, );

my $a= $ARGV[-1];
my $output= $1 if ( $a=~ /^([^\.]+)\./ );
`rm ./$output.HeatMap.svg` if (-e "$output.HeatMap.svg");

my ( $first_done, $lower ) = ( 0, 0 );

my ( @upper, @bottom );
my ( $distri_file_rev, $density_file_rev );
if ( @ARGV > 2 ){
	$density_file_rev= pop @ARGV;
	$distri_file_rev= pop @ARGV;
	$lower= 1;
	@upper= @ARGV;

        draw_distri_rev( $distri_file_rev );
        draw_density_rev( $density_file_rev );

	my $start_x= 0;
	while ( @upper > 0 ){
		my $tag= ( @upper == 2 ) ? 1 : 0;
		my $distri_file= shift @upper;
		my $density_file= shift @upper;

		if ( $first_done == 0 ){
			$start_x= draw_distri( $distri_file, 0, $tag );
			draw_density( $density_file, 0, $tag );
		}else{
                        draw_density( $density_file, $start_x, $tag );
			$start_x= draw_distri( $distri_file, $start_x, $tag );
		}
	}

}else{
	my $distri_file= shift;
        my $density_file= shift;
	
	draw_distri( $distri_file, 0, 1 );
	draw_density( $density_file, 0, 1 );
}

open OUT,">$output.HeatMap.svg" or die "$!";
print OUT $svg->xmlify();
#close OUT;

sub getMaxandAvr{
	my $Min_p = shift;
	my $Max_p = shift;
	my $Avr_p = shift;
	my $type_p = shift;
	my $type = shift;
	my $len = shift;
	my $total = 0;

	$Min_p->{$type} = 1;

	for (my $i=0; $i<$len; $i++)
	{
		if (($type_p->{$type}[$i] < $Min_p->{$type}) and ($type_p->{$type}[$i] > 0))
		{
			$Min_p->{$type} = $type_p->{$type}[$i];
		}
		elsif ($type_p->{$type}[$i] > $Max_p->{$type})
		{
			$Max_p->{$type} = $type_p->{$type}[$i];
		}

		$total += $type_p->{$type}[$i];
	}

#	$Min_p->{$type} = (int($Min_p->{$type}*100))/100;
#	$Max_p->{$type} = (int($Max_p->{$type}*100))/100;
	$Avr_p->{$type} = ((int($total/$len*100))/100 > 0) ? ((int($total/$len*100))/100) : 0.01;
	
}

sub draw_distri{
        my $file= shift;
	my $start_x= shift;
	my $tag= shift;
        my %color= ( "RTs"=>"#4E91C8", "DNA-TEs"=>"#FEEB84", "introns"=>"#A0CEE8", "exons"=>"#4FC3D0", );
        my %type;
	my ( $chr, $chr_len );
        open IN, $file or die "Can't open $file: $!\n";
        while(<IN>){
                chomp;
                my @t= split;
		$chr= $t[0];
                push @{$type{$t[1]}}, $t[3];
		$chr_len= $t[4];
        }
        close IN;

	my %Min;
	my %Max;
	my %Avr;
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "exons", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "introns", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "DNA-TEs", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "RTs", $chr_len);

	my ( $x, $y, $height ) = ( ($start_x+$edge_dis), $edge_dis, 100 );
	my $width= $chr_len*2;
	my $p= 0;
	$svg->rect('x',$x,'y',$y,'width',$width,'height',$height,'fill','#D3D3D3');
	for ( my $pos_x=$x; $pos_x<=($x+$width); $pos_x+=2 ){
		my $pos_y= $y + $height - $type{"exons"}->[$p];
		#my $height1= $type{"exons"}->[$p];
		my $height1= ($type{"exons"}->[$p]>1)?$type{"exons"}->[$p]:1;
		$svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$height1,'fill',$color{"exons"});
		$pos_y= $pos_y - $type{"introns"}->[$p];
		my $height2= ($type{"introns"}->[$p]>1)?$type{"introns"}->[$p]:1;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$height2,'fill',$color{"introns"});
		$pos_y= $pos_y - $type{"DNA-TEs"}->[$p];
		my $height3= ($type{"DNA-TEs"}->[$p]>1)?$type{"DNA-TEs"}->[$p]:1;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$height3,'fill',$color{"DNA-TEs"});
		$pos_y= $pos_y - $type{"RTs"}->[$p];
		my $height4= ($type{"RTs"}->[$p]>1)?$type{"RTs"}->[$p]:1;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$height4,'fill',$color{"RTs"});
		$p++;
	}
	my $chr_x= $start_x + $edge_dis + 3;
	my $chr_y= $edge_dis + 27;
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',24,'fill','black');
#	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',24,'stroke','black','fill','black');

	if ( $tag == 1 ){
		my ( $legend_w, $legend_h ) = ( 40, 15 );
		my $legend_x= $start_x + $edge_dis + $width + 10;
		my $legend_y= $edge_dis;
		my $text_x= $legend_x + $legend_w + 10;
		my $text_y= $legend_y + $legend_h;
		my $text_s= 20;
		my $data_x = $text_x + 120;
		my $data_y = $text_y;
		my $avr_x = $data_x + 100;
		my $avr_y = $text_y;

		my $line_x1 = $data_x - 5;
		my $line_y1 = $legend_y - 2;
		my $line_x2 = $avr_x + 50;
		my $line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);

		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"RTs"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"RTs",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"RTs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$Min{RTs} = (int($Min{RTs}*100))/100;
		$Max{RTs} = (int($Max{RTs}*100))/100;
#		$Avr{RTs} = (int($Avr{RTs}*100))/100;

		$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{RTs} - $Max{RTs}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{RTs} - $Max{RTs}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{RTs}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{RTs}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');

		$line_y1 = $text_y + 2;
		$line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1+97,'y1',$line_y1,'x2',$line_x1+97,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);

		$legend_y= $legend_y + $legend_h + 5;
		$text_y= $legend_y + $legend_h;
		$data_y = $text_y;
		$avr_y = $text_y;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"DNA-TEs"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"DNA-TEs",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
		$Min{"DNA-TEs"} = (int($Min{"DNA-TEs"}*100))/100;
		$Max{"DNA-TEs"} = (int($Max{"DNA-TEs"}*100))/100;
#		$Avr{"DNA-TEs"} = (int($Avr{"DNA-TEs"}*100))/100;
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"DNA-TEs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{\"DNA-TEs\"} - $Max{\"DNA-TEs\"}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{\"DNA-TEs\"} - $Max{\"DNA-TEs\"}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{\"DNA-TEs\"}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{\"DNA-TEs\"}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$line_y1 = $text_y + 2;
		$line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1+97,'y1',$line_y1,'x2',$line_x1+97,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
	
		$legend_y= $legend_y + $legend_h + 5;
		$text_y= $legend_y + $legend_h;
		$data_y = $text_y;
		$avr_y = $text_y;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"introns"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (introns)",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
		$Min{"introns"} = (int($Min{"introns"}*100))/100;
		$Max{"introns"} = (int($Max{"introns"}*100))/100;
#		$Avr{"introns"} = (int($Avr{"introns"}*100))/100;
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (introns)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{introns} - $Max{introns}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{introns} - $Max{introns}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{introns}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{introns}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$line_y1 = $text_y + 2;
		$line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1+97,'y1',$line_y1,'x2',$line_x1+97,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
	
		$legend_y= $legend_y + $legend_h + 5;
		$text_y= $legend_y + $legend_h;
		$data_y = $text_y;
		$avr_y = $text_y;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"exons"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (exons)",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
		$Min{"exons"} = (int($Min{"exons"}*100))/100;
		$Max{"exons"} = (int($Max{"exons"}*100))/100;
#		$Avr{"exons"} = (int($Avr{"exons"}*100))/100;
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (exons)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{exons} - $Max{exons}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{exons} - $Max{exons}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{exons}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
	#	$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{exons}",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
		$line_y1 = $text_y + 2;
		$line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1+97,'y1',$line_y1,'x2',$line_x1+97,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
		
		$data_y = $data_y + 20;
		$data_x += 20;
		$svg->text('x',$data_x,'y',$data_y,'-cdata',"range",'font-family',"ArialNarrow",'font-size',18,'fill','black');
		$avr_x += 10;
		$svg->text('x',$avr_x,'y',$data_y,'-cdata',"avg",'font-family',"ArialNarrow",'font-size',18,'fill','black');
#		$line_y1 = $text_y + 2;
#               $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);

		$line_y1 = $data_y + 5;
		$line_y2 = $line_y1;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
                $svg->line('x1',$line_x1+97,'y1',$line_y1,'x2',$line_x1+97,'y2',$legend_y - 2,'stroke','grey','stroke-width',1);
		$line_x1 = $legend_x;
                $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','grey','stroke-width',1);
	}
	
	$first_done= 1;
	$start_x= $start_x + $width + 20;
	return ( $start_x );
}


sub draw_density{
	my $file= shift;
	my $start_x= shift;
	my $tag= shift;
	my %type;
	my $chr;
	my $chr_len;
	#my %max;
	
	my ( @color_arr, @arr );
	@arr= (0,0,128);
	for ( my $i=128; $i<255; $i++ ){
		push @color_arr, @arr;
		$arr[2]++;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]++;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[0]++;
		$arr[2]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<127; $i++ ){
		$arr[0]--;
		push @color_arr, @arr;
	}
	open IN, $file or die "Can't open $file: $!\n";
	while(<IN>){
		chomp;
		my @t= split /\t/;
		$chr= $t[0];
		push @{$type{$t[1]}}, $t[3];
		$chr_len= $t[4];
	}
	
	my $pos_y= $edge_dis + 106;
	my $type_h= 15;
	my $text_s = 20;
	my ( $text_x, $text_y, $text_x2, $text_y2);

#	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)", "Paralogs (exons)" );
	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)");
	my %nickName = ("LTR-RT/gypsy"=>"LTR-RT/ gypsy", "LTR-RT/copia"=>"LTR-RT/ copia", "DNA-TE/CACTA"=>"DNA-TE/ CACTA",
				"DNA-TE/MITE"=>"DNA-TE/ MITE", "Genes (exons)"=>"Genes");
	my %Min;
	my %Max;
	my %Avr;
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "LTR-RT/gypsy", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "LTR-RT/copia", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "DNA-TE/CACTA", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "DNA-TE/MITE", $chr_len);
	getMaxandAvr(\%Min, \%Max, \%Avr, \%type, "Genes (exons)", $chr_len);

	my ($data_x, $data_y, $avr_x, $avr_y);

	foreach my $tp ( @type_arr ){
		my $pos_x= $start_x + $edge_dis;
#print "$tp: $Max{$tp}\n";
		for ( my $i=0; $i<@{$type{$tp}}; $i++ ){
		#	my $j= int($type{$tp}->[$i]*10.18)*3;
			my $j= int($type{$tp}->[$i]*100*10.18/$Max{$tp})*3;
			my $r_color= $color_arr[$j];
			my $g_color= $color_arr[$j+1];
			my $b_color= $color_arr[$j+2];
			$svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$type_h,'fill',"rgb($r_color,$g_color,$b_color)");
			$pos_x+= 2;
		}
		$text_x= $pos_x + 10;
		$text_x2 = $pos_x + 330;
		$text_y= $pos_y + $text_s - 5;
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'stroke','black','fill','black') if ( $tag == 1 );
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"$nickName{$tp}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black') if ( $tag == 1 );
		$data_x = $text_x + 170;
		$data_y = $text_y;
		$avr_x = $data_x + 100;
		$avr_y = $data_y;

#		$Min{$tp} = ((int($Min{$tp}*100))/100 > 0)?((int($Min{$tp}*100))/100):0.01;
#		$Max{$tp} = ((int($Max{$tp}*100))/100 > 0)?((int($Max{$tp}*100))/100):0.01;
		$Min{$tp} = (int($Min{$tp}*100))/100;
		$Max{$tp} = (int($Max{$tp}*100))/100;
		$svg->text('x',$data_x,'y',$data_y,'-cdata',"$Min{$tp} - $Max{$tp}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
		$svg->text('x',$avr_x,'y',$avr_y,'-cdata',"$Avr{$tp}",'font-family',"ArialNarrow",'font-size',$text_s,'fill','black');
		

		$pos_y= $pos_y + $type_h + 5;

                $svg->line('x1',$text_x,'y1',$text_y+2,'x2',$text_x2,'y2',$text_y+2,'stroke','black','stroke-width',1);
	}

	my $line_x1 = $start_x + $edge_dis + $chr_len*2 + 175;
	my $line_x2 = $line_x1 + 97;
	my $line_x3 = $line_x1 + 155;
	my $line_y1 = $edge_dis + 100;
	my $line_y2 = $pos_y - 3;

#print "x1: $line_x1\tx2: $line_x2\tx3: $line_x3\n";
        $svg->line('x1',$line_x1,'y1',$line_y1,'x2',$line_x1,'y2',$line_y2,'stroke','black','stroke-width',1);
        $svg->line('x1',$line_x2,'y1',$line_y1,'x2',$line_x2,'y2',$line_y2,'stroke','black','stroke-width',1);
        $svg->line('x1',$line_x3,'y1',$line_y1,'x2',$line_x3,'y2',$line_y2,'stroke','black','stroke-width',1);

	for ( my $i=$edge_dis; $i<=($edge_dis+2*$chr_len); $i+=200 ){
                my $length= ( $i-$edge_dis )/20;
                my $x1= ( $i==$edge_dis ) ? $i : ( $i-18 );
		my $y1= $text_y + 25;
                $svg->text('x',$x1,'y',$y1,'-cdata',$length."M",'font-family',"ArialNarrow",'font-size',20,'fill','black');
               # $svg->text('x',$x1,'y',$y1,'-cdata',$length."M",'font-family',"ArialNarrow",'font-size',20,'stroke','black','fill','black');
                $svg->line('x1',$i,'y1',$text_y,'x2',$i,'y2',$text_y+4,'stroke','black','stroke-width',2);
#print "i:$i; text_y: $text_y; length:$length; chr len:$chr_len\n";
        }


	if ( $lower == 0 ){
		my $legend_x= $text_x + 340;
		my $legend_y= $text_y;
		my $legend_w= 16;
		for ( my $i=0; $i<=($text_y-$edge_dis); $i+=1 ){
			my $j= int($i*1019/($text_y-$edge_dis+1))*3;
			my $r_color= $color_arr[$j];
			my $g_color= $color_arr[$j+1];
			my $b_color= $color_arr[$j+2];
			$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',1,'fill',"rgb($r_color,$g_color,$b_color)");
			$legend_y= $text_y - $i;
		}

		$text_x= $legend_x + 35;
		$text_y= $text_y - ($text_y-$edge_dis)/5;
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',20,'fill','black','transform',"rotate(-90,$text_x,$text_y)");
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',20,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
	}
	
}


=pod sub draw_distri_rev{
        my $file= shift;
        my %color= ( "RTs"=>"#4E91C8", "DNA-TEs"=>"#FEEB84", "introns"=>"#A0CEE8", "exons"=>"#4FC3D0", );
        my %type;
	my ( $chr, $chr_len );
        open IN, $file or die "Can't open $file: $!\n";
        while(<IN>){
                chomp;
                my @t= split;
		$chr= $t[0];
                push @{$type{$t[1]}}, $t[3];
		$chr_len= $t[4];
        }
        close IN;
	my ( $x, $y, $height ) = ( $edge_dis, 306, 50 );
	my $width= $chr_len;
	my $p= 0;
	$svg->rect('x',$x,'y',$y,'width',$width,'height',$height,'fill','#D3D3D3');
	for ( my $pos_x=$x; $pos_x<=($x+$width); $pos_x+=1 ){
		my $pos_y= $y + $height - $type{"exons"}->[$p]/2;
		my $height1= $type{"exons"}->[$p]/2;
		$svg->rect('x',$pos_x,'y',$pos_y,'width',1,'height',$height1,'fill',$color{"exons"});
		$pos_y= $pos_y - $type{"introns"}->[$p]/2;
		my $height2= $type{"introns"}->[$p]/2;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',1,'height',$height2,'fill',$color{"introns"});
		$pos_y= $pos_y - $type{"DNA-TEs"}->[$p]/2;
		my $height3= $type{"DNA-TEs"}->[$p]/2;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',1,'height',$height3,'fill',$color{"DNA-TEs"});
		$pos_y= $pos_y - $type{"RTs"}->[$p]/2;
		my $height4= $type{"RTs"}->[$p]/2;
                $svg->rect('x',$pos_x,'y',$pos_y,'width',1,'height',$height4,'fill',$color{"RTs"});
		$p++;
	}
	my $chr_x= $edge_dis + 3;
	my $chr_y= 306 + 15;
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',12,'stroke','black','fill','black');

	my ( $legend_w, $legend_h ) = ( 15, 5 );
	my $legend_x= $edge_dis + $width + 5;
	my $legend_y= 306;
	my $text_x= $legend_x + $legend_w + 5;
	my $text_y= $legend_y + $legend_h + 2;
	my $text_s= 8;

	$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"RTs"});
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"RTs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');

	$legend_y= $legend_y + 2*$legend_h;
	$text_y= $legend_y + $legend_h + 2;
	$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"DNA-TEs"});
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"DNA-TEs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
	
	$legend_y= $legend_y + 2*$legend_h;
	$text_y= $legend_y + $legend_h + 2;
	$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"introns"});
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (introns)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
	
	$legend_y= $legend_y + 2*$legend_h;
	$text_y= $legend_y + $legend_h + 2;
	$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"exons"});
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (exons)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
	
#	open OUT,">>$chr.HeatMap.svg" or die "$!";
#	print OUT $svg->xmlify();
	#close OUT;
}


sub draw_density_rev{
	my $file= shift;
	my %type;
	my $chr;
	#my %max;
	
	my ( @color_arr, @arr );
	@arr= (0,0,128);
	for ( my $i=128; $i<255; $i++ ){
		push @color_arr, @arr;
		$arr[2]++;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]++;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[0]++;
		$arr[2]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<255; $i++ ){
		$arr[1]--;
		push @color_arr, @arr;
	}
	for ( my $i=0; $i<127; $i++ ){
		$arr[0]--;
		push @color_arr, @arr;
	}
	open IN, $file or die "Can't open $file: $!\n";
	while(<IN>){
		chomp;
		my @t= split /\t/;
		$chr= $t[0];
		push @{$type{$t[1]}}, $t[3];
	}
	
	my $pos_y= 228;
	my $type_h= 10;
	my ( $text_x, $text_y );
	#my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)", "Paralogs (exons)" );
#	my @type_arr= ( "Paralogs (exons)", "Genes (exons)", "DNA-TE/MITE", "DNA-TE/CACTA", "LTR-RT/copia", "LTR-RT/gypsy" );
	my @type_arr= ( "genes (exons)", "DNA-TE/MITE", "DNA-TE/CACTA", "LTR-RT/copia", "LTR-RT/gypsy" );
	foreach my $tp ( @type_arr ){
		my $pos_x= $edge_dis;
		for ( my $i=0; $i<@{$type{$tp}}; $i++ ){
			my $j= int($type{$tp}->[$i]*10.18)*3;
			my $r_color= $color_arr[$j];
			my $g_color= $color_arr[$j+1];
			my $b_color= $color_arr[$j+2];
			$svg->rect('x',$pos_x,'y',$pos_y,'width',1,'height',$type_h,'fill',"rgb($r_color,$g_color,$b_color)");
			$pos_x+= 1;
		}
		$text_x= $pos_x + 5;
		$text_y= $pos_y + $type_h;
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'fill','black');
	#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'stroke','black','fill','black');
		$pos_y= $pos_y + $type_h + 3;
	}

	my $legend_x= $text_x + 120;
	my $legend_y= 356;
	my $legend_w= 8;
	for ( my $i=0; $i<=(356-228); $i+=0.2 ){
		my $j= int($i*1019/(356-228+1))*3;
		my $r_color= $color_arr[$j];
		my $g_color= $color_arr[$j+1];
		my $b_color= $color_arr[$j+2];
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',0.2,'fill',"rgb($r_color,$g_color,$b_color)");
		$legend_y= 356 - $i;
	}

	$text_x= $legend_x + 18;
	$text_y= $text_y - ($text_y-$edge_dis)/4;
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',8,'fill','black','transform',"rotate(-90,$text_x,$text_y)");
#	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',8,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
	
	
}
=cut
