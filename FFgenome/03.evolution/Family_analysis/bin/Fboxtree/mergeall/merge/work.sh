perl merge.pl --os ../rice/4tree --ob ../brachyantha/4tree --sb ../sorghum/4tree
perl cleanUnknown.pl --dir ./
perl checkoverlap.pl --dir ./ > overlap
perl assignPfam.pl --dir ./

