use Bio::MapIO;
use Bio::Map::Physical;
use Bio::Map::Clone;
use Bio::Map::Contig;
use Bio::Map::FPCMarker;
use Bio::Range;

my $mapio=Bio::MapIO->new(-format => "fpc",-file => "OB__B.fpc",
                               -readcor => 0);

my $physical = $mapio->next_map();
foreach my $clone ($physical->each_cloneid()){
    #print "$clone\n";
    my $cloneobj=$physical->get_cloneobj($clone);
    my $name=$cloneobj->name();
    my $contig=$cloneobj->contigid();
    print "$clone\t$name\t$contig\n";   
}
