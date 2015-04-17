## cut a long sequence into 40 fregment with length between 200K-2M.


while ($i<40){
srand;
my $int1=int(rand(3000));
my $int2=int(rand(900));

my $start=$int1*10000+1; ##start point on chromosome with length 30M
my $length=($int2+100)*2000; ##length of sequence from 200K-2M
my $end=$start+$length;
system "perl getsegment.pl $start $end chr04.con";
$i++;
}
system "cat *.txt > psudoseq";
system "rm *.txt";
