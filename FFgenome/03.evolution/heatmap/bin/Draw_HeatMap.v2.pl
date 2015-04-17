#!/usr/bin/perl
#author: Li zhenyu (lizhenyu@genomics.cn)
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
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',24,'stroke','black','fill','black');

	if ( $tag == 1 ){
		my ( $legend_w, $legend_h ) = ( 45, 15 );
		my $legend_x= $start_x + $edge_dis + $width + 10;
		my $legend_y= $edge_dis;
		my $text_x= $legend_x + $legend_w + 10;
		my $text_y= $legend_y + $legend_h + 5;
		my $text_s= 20;

		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"RTs"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"RTs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');

		$legend_y= $legend_y + $legend_h + 10;
		$text_y= $legend_y + $legend_h + 2;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"DNA-TEs"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"DNA-TEs",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
	
		$legend_y= $legend_y + $legend_h + 10;
		$text_y= $legend_y + $legend_h + 2;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"introns"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (introns)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
	
		$legend_y= $legend_y + $legend_h + 10;
		$text_y= $legend_y + $legend_h + 2;
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',$legend_h,'fill',$color{"exons"});
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"gene (exons)",'font-family',"ArialNarrow",'font-size',$text_s,'stroke','black','fill','black');
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
#	A01     LTR-RT/gypsy    0       1.4268  258
	open IN, $file or die "Can't open $file: $!\n";
	while(<IN>){
		chomp;
		my @t= split /\t/;
		$chr= $t[0];
		push @{$type{$t[1]}}, $t[3];
		$chr_len= $t[4];
	}
	
	my $pos_y= $edge_dis + 106;
	my $type_h= 20;
	my ( $text_x, $text_y );
#	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)", "Paralogs (exons)" );
	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)");
	foreach my $tp ( @type_arr ){
		my $pos_x= $start_x + $edge_dis;
		for ( my $i=0; $i<@{$type{$tp}}; $i++ ){
			my $j= int($type{$tp}->[$i]*22.18)*3;
			my $r_color= $color_arr[$j];
			my $g_color= $color_arr[$j+1];
			my $b_color= $color_arr[$j+2];
			$svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$type_h,'fill',"rgb($r_color,$g_color,$b_color)");
			$pos_x+= 2;
		}
		$text_x= $pos_x + 10;
		$text_y= $pos_y + $type_h;
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'stroke','black','fill','black') if ( $tag == 1 );
		$pos_y= $pos_y + $type_h + 5;
	}

	for ( my $i=$edge_dis; $i<=($edge_dis+2*$chr_len); $i+=200 ){
                my $length= ( $i-$edge_dis )/20;
                my $x1= ( $i==$edge_dis ) ? $i : ( $i-18 );
		my $y1= $text_y + 25;
                $svg->text('x',$x1,'y',$y1,'-cdata',$length."M",'font-family',"ArialNarrow",'font-size',20,'stroke','black','fill','black');
                $svg->line('x1',$i,'y1',$text_y,'x2',$i,'y2',$text_y+4,'stroke','black','stroke-width',2);
#print "i:$i; text_y: $text_y; length:$length; chr len:$chr_len\n";
        }


	if ( $lower == 0 ){
		my $legend_x= $text_x + 200;
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
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',20,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
	}
	
}


sub draw_distri_rev{
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
	my @type_arr= ( "Paralogs (exons)", "Genes (exons)", "DNA-TE/MITE", "DNA-TE/CACTA", "LTR-RT/copia", "LTR-RT/gypsy" );
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
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',8,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
	
	
}
=cut
