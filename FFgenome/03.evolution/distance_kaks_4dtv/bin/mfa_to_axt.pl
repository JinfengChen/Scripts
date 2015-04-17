#!/usr/bin/perl
use strict;

##Author: Fan Wei
##Email: fanw@genomic.org.cn
##Data: 2008-12-22

my $file = shift;

mfa2axt($file);

##get matrix of Ka and Ks from KaKs_caculator's result
sub mfa2axt{
        my $mfa_file = shift;
        my (%name_seq,%pair,$output);

        open OUT, ">$mfa_file.axt" || die "fail $mfa_file.axt";

        open IN, $mfa_file || die "fail open $mfa_file\n";
        $/=">"; <IN>; $/="\n";
        while (<IN>) {
                my $name = $1 if(/^(\S+)/);
                $/=">";
                my $seq = <IN>;
                chomp $seq;
                $seq =~ s/\s//g;
                $/="\n";
                $name_seq{$name} = $seq;
        }
        close IN;

        foreach my $first (sort keys %name_seq) {
                foreach my $second (sort keys %name_seq) {
                        next if($first eq $second || exists $pair{"$second&$first"});
                        $pair{"$first&$second"} = 1;
                }
        }

        foreach (sort keys %pair) {
                if (/([^&]+)&([^&]+)/) {
                        print OUT $_."\n".$name_seq{$1}."\n".$name_seq{$2}."\n\n";
                }
        }

        close OUT;
}
