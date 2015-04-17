#!/bin/bash -eu

echo "*** Performing GFF Merging ***"

if [ $# -lt 2 ] ## if argv < 2 
then
	echo "Usage: $0 <output> <GFF file to merge>..."
	exit 1
fi

output=`cd \`dirname $1\`; pwd`/`basename $1`  ##output file, which is first argv, $1
tempout=$output.tmp
input=${output/.gff/}.raw.gff                  ##temp input file


shift


## do argv input file and put all content in temp input file
inputs=''
for i in $* 
do
	inputs="$inputs `cd \`dirname $i\`; pwd`/`basename $i`"
done

cat $inputs > $input

## collect features and sources of SV callers
features="`cut -f 3 $input | sort -u`"
sources="`cut -f 2 $input | sort -u`"


## do merge of gff
for feat in $features
do
	f=$input.${feat}
	grep "	$feat	" $input  > $f  ## first, split each features into one file, test.raw.gff.Deletion
	for src in $sources 
	do
		grep "$src" $f | mergeBed -i stdin | awk '{print $1"\t"$2"\t"$3"\t""'$src'"}' > $f.$src ## second, split each source of each features and merge, test.raw.gff.Deletion.BD
                
	done
	cat $f.* | msort -k 1,n2 > $f   ## third, combine source of each features
	> $f.dup
	> $f.uni
	intersectBed -a $f -b $f -f 0.5 -r -c | awk '{print $0 > ($NF>1? "'$f.dup'": "'$f.uni'")}' ## do overlap
	mergeBed -n -nms -i $f.uni > $f.uni.merged ## merge unique
	mergeBed -n -nms -i $f.dup > $f.dup.merged ## merge overlap
	for i in $f.*.merged ## parse merged file into final gff 
	do
		dup="0"
		if [ -n "${i/*.uni.merged/}" ]; then dup="1"; fi ## have overlap, so should be high confidential
		awk -F '\t' '{qual="LowQual"; if ('$dup' && $4>=2) qual="PASS"; print $1"\tHugeSeq\t'$feat'\t"$2"\t"$3"\t"qual"\t.\t.\tEVENTS="$5";Caller="$4}' $i
	done
	rm $f $f.*
done > $tempout

msort -k 1,n4 $tempout > $output
rm $tempout $input 

