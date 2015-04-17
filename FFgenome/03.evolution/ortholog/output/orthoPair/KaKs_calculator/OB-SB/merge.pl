chomp $ARGV[0];
my $dir=$ARGV[0];
my @inf=glob("$dir/*.inf");
foreach(@inf){
    `cat $_ >> $dir/KaKs.summary`;
    `rm $_`;
}
