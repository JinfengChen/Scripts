#!/usr/bin/perl
#author: Li zhenyu (lizhenyu@genomics.cn)
# usage: 
# 	(1) draw each chromosome separately: perl $0 -dir data_directory &
#	(2) draw correlated chromosome: perl $0 -dir data_directory -pos all.gene.postion -block connect.txt
#use strict;
use Data::Dumper;
use lib "$Bin/../lib";
use SVG;

my ( $newwidth, $newheight, $edge_dis ) = ( 1200, 800, 50 );
my $svg=SVG->new( 'width',$newwidth,'height',$newheight, );

my %related;
my $related_file= pop @ARGV;
read_related( $related_file, \%related ) if ( $related_file ne "" );

my %line;
my $line_file= pop @ARGV;
read_line( $line_file, \%line );

my $a= $ARGV[-2];
my $output= $1 if ( $a=~ /^([^\.]+)\./ );
`rm ./$output.HeatMap.svg`;

my ( $first_done, $lower ) = ( 0, 0 );

my ( @upper, @bottom );
my ( $related_s, $related_e, $target_s, $target_e );
my ( $distri_file_rev, $density_file_rev );
if ( @ARGV > 2 ){
	$density_file_rev= pop @ARGV;
	$distri_file_rev= pop @ARGV;
	$lower= 1;
	@upper= @ARGV;

        draw_distri_rev( $distri_file_rev );
        draw_density_rev( $density_file_rev );

	my $start_x= 0;
	my $i= 0;
	while ( @upper > 0 ){
		my $tag= ( @upper == 2 ) ? 1 : 0;
		my $distri_file= shift @upper;
		my $density_file= shift @upper;
		
		my $related_chr= $1 if ( $distri_file=~ /^([^\.]+)\./ );
		$related_s= $related{$related_chr}->[2];
		$related_e= $related{$related_chr}->[3];
		$target_s= $related{$related_chr}->[0];
		$target_e= $related{$related_chr}->[1];
		
#print "related chr:$related_chr; related s:$related_s; related e:$related_e; distri_file: $distri_file; density_file:$density_file\n";
		if ( $first_done == 0 ){
			$start_x= draw_distri( $distri_file, 0, $tag, $related_s, $related_e, $target_s, $target_e );
			draw_density( $density_file, 0, $tag, $related_s, $related_e, $target_s, $target_e  );
		}else{
                       	draw_density( $density_file, $start_x, $tag, $related_s, $related_e, $target_s, $target_e );
			$start_x= draw_distri( $distri_file, $start_x, $tag, $related_s, $related_e, $target_s, $target_e );
		}
		$i++;
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
	my $related_start= shift;
	my $related_end= shift;
	my $target_start= shift;
	my $target_end= shift;
	( $related_start, $related_end ) = ( $related_start/100000, $related_end/100000 );
	( $target_start, $target_end ) = ( $target_start/100000, $target_end/100000 );
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
	#my ( $x, $y, $height ) = ( ($start_x+$edge_dis+$target_start), $edge_dis, 50 );
	#my ( $x, $y, $height ) = ( ($start_x+$edge_dis+$related_start), $edge_dis, 50 );
	my ( $x, $y, $height ) = ( ($start_x+$edge_dis), $edge_dis, 100 );
	#my $width= $chr_len;
	my $width= ( ($related_end - $related_start) >0 ) ? 2*($related_end - $related_start) : $chr_len*2;
	my $p= ( int($related_start)>0 ) ? int($related_start): 0;
	$svg->rect('x',$x,'y',$y,'width',$width,'height',$height,'fill','#D3D3D3');
print "related_start:$related_start;related_end:$related_end;width:$width;p:$p; x:$x\n";
	for ( my $pos_x=$x; $pos_x<=($x+$width); $pos_x+=2 ){
		my $pos_y= $y + $height - $type{"exons"}->[$p];
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

	my $scalar;
=pod	if ( $width > 200 ){
		$scalar= 100;
	}elsif( $width > 50 ){
		$scalar= 50;
	}else{
		$scalar= 30;
	}
=cut
	for ( my $i=$x; $i<=($x+$width); $i+= $width){
		my $length= ( $i==$x ) ? (int($related_start ))/10 : (int($related_end))/10;
		my $x1= ( $i==$x )? $i : ( $i-45 );
print "length:$length\n";
		$svg->text('x',$x1,'y',46,'-cdata',$length."M",'font-family',"ArialNarrow",'font-size',18,'stroke','black','fill','black');
		$svg->line('x1',$i,'y1',$edge_dis,'x2',$i,'y2',$edge_dis-4,'stroke','black','stroke-width',2);
	}

	my $chr_x= $x + 3;
	my $chr_y= $edge_dis + 27;
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',24,'stroke','black','fill','black');

	if ( $tag == 1 ){
		my ( $legend_w, $legend_h ) = ( 45, 15 );
		my $legend_x= $start_x + $edge_dis + $width + 10;
		#my $legend_x= $start_x + $edge_dis + $target_start + $related_end - $related_start + 5;
		my $legend_y= $edge_dis;
		my $text_x= $legend_x + $legend_w + 5;
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
	#$start_x= ( $related_start > 0 ) ? 0 : ($start_x + $width + 20 );
	$start_x= ( $related_start > 0 ) ? ( $start_x + $width + 5 ) : ($start_x + $width + 20 );
#print "stsrt_x:$start_x\n";
	return ( $start_x );
}


sub draw_density{
	my $file= shift;
	my $start_x= shift;
	my $tag= shift;
	my $related_start= shift;
	my $related_end= shift;
	my $target_start= shift;
        my $target_end= shift;
        ( $related_start, $related_end ) = ( $related_start/100000, $related_end/100000 );
        ( $target_start, $target_end ) = ( $target_start/100000, $target_end/100000 );

	my %type;
	my ( $chr, $chr_len );
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
	my $type_h= 20;
	my ( $text_x, $text_y );
	my @type_arr= ( "LTR-RT/gypsy", "LTR-RT/copia", "DNA-TE/CACTA", "DNA-TE/MITE", "Genes (exons)", "Paralogs (exons)" );
	my $pos_x;
	my $x;
	foreach my $tp ( @type_arr ){
		my $p= ( $related_start > 0 ) ? int($related_start) : 0;
		my $width= ( ($related_end - $related_start) >0 ) ? 2*($related_end - $related_start) : $chr_len*2;
		my $pos_end= ( $related_end > 0 ) ? int($related_end) : (@{$type{$tp}}-1);
		#my $x= $start_x + $edge_dis + int($target_start);
		$x= $start_x + $edge_dis;
		$pos_x= $x;
		#my $x= $start_x + $edge_dis + $pos_start;
		for ( $i=$p; $i<=$pos_end; $i++ ){
		#for ( $pos_x=$x; $pos_x<=($x+$width); $pos_x+=2 ){
			my $j= int($type{$tp}->[$i]*10.18)*3;
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

	for ( my $i=0; $i<@{$line{$chr}}; $i++ ){
		my $ln_p_x1= $edge_dis + 2*$line{$chr}->[$i][0]/100000;
		my $ln_p_x2= $x + 2*($line{$chr}->[$i][1]/100000 - $related_start);
		my $ln_p_y1= $text_y + 50 - 1;
		my $ln_p_y2= $text_y + 1;
		$svg->line('x1',$ln_p_x1,'y1',$ln_p_y1,'x2',$ln_p_x2,'y2',$ln_p_y2,'stroke','#D3D3D3','stroke-width',1);
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
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',8,'stroke','black','fill','black','transform',"rotate(-90,$text_x,$text_y)");
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
	my ( $x, $y, $height ) = ( $edge_dis, 502 , 100 );
	my $width= $chr_len*2;
	my $p= 0;
	$svg->rect('x',$x,'y',$y,'width',$width,'height',$height,'fill','#D3D3D3');
	for ( my $pos_x=$x; $pos_x<=($x+$width); $pos_x+=2 ){
		my $pos_y= $y + $height - $type{"exons"}->[$p];
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
	my $chr_x= $edge_dis + 3;
	my $chr_y= 502 + 27;
	$svg->text('x',$chr_x,'y',$chr_y,'-cdata',"$chr",'font-family',"ArialNarrow",'font-size',24,'stroke','black','fill','black');

	my ( $legend_w, $legend_h ) = ( 45, 15 );
	my $legend_x= $edge_dis + $width + 10;
	my $legend_y= 502;
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
	
	for ( my $i=$edge_dis; $i<=($edge_dis+$width); $i+=200 ){
		my $length= ( $i-50 )/20;
		my $x1= ( $i==$edge_dis )? $i : ( $i-14 );
		$svg->text('x',$x1,'y',627,'-cdata',$length."M",'font-family',"ArialNarrow",'font-size',18,'stroke','black','fill','black');
		$svg->line('x1',$i,'y1',606,'x2',$i,'y2',602,'stroke','black','stroke-width',2);
	}
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
	
	my $pos_y= 351;
	my $type_h= 20;
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
			$svg->rect('x',$pos_x,'y',$pos_y,'width',2,'height',$type_h,'fill',"rgb($r_color,$g_color,$b_color)");
			$pos_x+= 2;
		}
		$text_x= $pos_x + 5;
		$text_y= $pos_y + $type_h;
		$svg->text('x',$text_x,'y',$text_y,'-cdata',"$tp",'font-family',"ArialNarrow",'font-size',$type_h,'stroke','black','fill','black');
		$pos_y= $pos_y + $type_h + 5;
	}

	my $legend_x= $text_x + 250;
	my $legend_y= 602;
	my $legend_w= 16;
	for ( my $i=0; $i<=(602-351); $i+=1 ){
		my $j= int($i*1019/(602-351+1))*3;
		my $r_color= $color_arr[$j];
		my $g_color= $color_arr[$j+1];
		my $b_color= $color_arr[$j+2];
		$svg->rect('x',$legend_x,'y',$legend_y,'width',$legend_w,'height',1,'fill',"rgb($r_color,$g_color,$b_color)");
		$legend_y= 602 - $i;
	}

	$text_x= $legend_x + 35;
	$text_y= $text_y + 50;
	#$text_y= $text_y - ($text_y-$pos_y)/2;
	$svg->text('x',$text_x,'y',$text_y,'-cdata',"% bp [0 to max]",'font-family',"ArialNarrow",'font-size',20,'stroke','black','fill','black','transform',"rotate(270,$text_x,$text_y)");
	
	
}


sub read_related{
	my $file= shift;
	my $related_p= shift;

	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
		chomp;
		my @t= split;
		push @{$related_p->{$t[3]}}, ($t[1], $t[2], $t[4], $t[5]);
#print "$t[3]\n";
	}
}

sub read_line{
	my $file= shift;
	my $line_p= shift;

	open IN, $file or die "Can't open $file:$!\n";
	while(<IN>){
		chomp;
		my @t= split;
		my $target= $1 if ( $t[0] =~ /^[^\_]+\_\_([^\_]+)\_/ );
		push @{$line_p->{$target}}, [$t[2], $t[5]] if ( ($t[5]>=$related{$target}[2]) && ($t[5]<=$related{$target}[3]) );
	}
}
