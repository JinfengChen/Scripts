grep "Os01" colinear.lst > chr01.list
grep "Os02" colinear.lst > chr02.list
grep "Os03" colinear.lst > chr03.list
grep "Os04" colinear.lst > chr04.list
grep "Os05" colinear.lst > chr05.list
grep "Os06" colinear.lst > chr06.list
grep "Os07" colinear.lst > chr07.list
grep "Os08" colinear.lst > chr08.list
grep "Os09" colinear.lst > chr09.list
grep "Os10" colinear.lst > chr10.list
grep "Os11" colinear.lst > chr11.list
grep "Os12" colinear.lst > chr12.list

perl besthit.pl --table 3wayv1.4.blast &


