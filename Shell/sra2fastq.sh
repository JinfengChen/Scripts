#!/bin/sh

sra2fastq () {
    for file in `ls *.sra 2> md5.err`;do
       ((convertfile++))
       echo "Converting $file"
       /opt/sratoolkit/2.1.16/bin/fastq-dump --outdir $1 --split-files $file 
    done
    rm md5.err
    return 1
}


do_recursive() {   # recursive execute "checkmd5" in this directory

    sra2fastq $1
    for filename in `ls`;do
         if [ -d "$filename" ]
         then
             cd $filename
             do_recursive $1
             cd ..
         fi
    done
    return 1
}



#make new dir fastq and specify the fastq directory to this dir
root_dir=`pwd`
if [ -e "fastq" ]
then
   continue
else
   echo "Make new directory fastq"
   mkdir "fastq"
fi

target="$root_dir/fastq"
echo "Target directory of fastq file: $target"

#do converting
do_recursive $target
echo "Totally converted file: $convertfile"


