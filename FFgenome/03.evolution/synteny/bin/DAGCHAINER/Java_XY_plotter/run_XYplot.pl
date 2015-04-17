#!/usr/bin/env perl

use strict;

my $libPath = $0;
$libPath =~ s/[^\/]+$//;
$libPath .= "lib";

print STDERR "libpath: $libPath\n";

my $cmd = "java -Xmx300M -cp  $libPath  JsyntenyView @ARGV";

exec $cmd;

die "ERROR, couldn't exec $cmd\n";


 




