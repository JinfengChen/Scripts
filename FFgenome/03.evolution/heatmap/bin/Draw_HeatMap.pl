#!/usr/bin/perl
#use strict;
use Data::Dumper;
use lib "$Bin/../lib";
use SVG;

my $distri_file= shift;
my $density_file= shift;
my $output= $1 if ( $distri_file=~ /^([^\.]+)\./ );
my ( $newwidth, $newheight, $edge_dis ) = ( 800, 800, 50 );
my $svg=SVG->new( 'width',$newwidth,'height',$newheight, );

`rm ./$output.HeatMap.svg`;

draw_distri( $distri_file );
draw_density( $density_file );

open OUT,">$output.HeatMap.svg" or die "$!";
print OUT $svg->xmlify();
#close OUT;

sub draw_distri{
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
	my ( $x, $y, $height ) = ( $edge_dis, $edge_dis, 50 );
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
	my $chr_y= $edge_dis + 15;
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',12,'stroke','black','fill','black');

	my ( $legend_w, $legend_h ) = ( 15, 5 );
	my $legend_x= $edge_dis + $width + 5;
	my $legend_y= $edge_dis;
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


sub draw_density{
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
	
	my $pos_y= $edge_dis + 53;
	my $type_h= 10;
	my ( $text_x, $text_y );
	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)", "Paralogues" );
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
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'stroke','black','fill','black');
		$pos_y= $pos_y + $type_h + 3;
	}

	my $legend_x= $text_x + 100;
	my $legend_y= $text_y;
	my $legend_w= 8;
	for ( my $i=0; $i<=($text_y-$edge_dis); $i+=0.2 ){
		my $j= int($i*1019/($text_y-$edge_dis+1))*3;
		my $r_color= $color_arr[$j];
		my $g_color= $color_arr[$j+1];
		my $b_color= $color_arr[$j+2];
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',0.2,'fill',"rgb($r_color,$g_color,$b_color)");
		$legend_y= $text_y - $i;
	}

	$text_x= $legend_x + 18;
	$text_y= $text_y - ($text_y-$edge_dis)/4;
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',8,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
	
	
#LTR-RT/gypsy
#LTR-RT/copia
#DNA-TE/CACTA
#DNA-TE/MITE
#Genes (exons)
#Paralogues
}
