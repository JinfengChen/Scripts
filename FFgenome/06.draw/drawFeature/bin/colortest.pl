my @refcol=colortab();
my $number=@refcol;
foreach(@refcol){
  print "$_\n";
}
#print "$refcol[$number-1]\n";
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

