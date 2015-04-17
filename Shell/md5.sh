#!/bin/sh
#make/clean/check md5 code for file in input directory
#Usege: ./md5.sh make|clean|check suffix


makemd5 () { # make md5 for *.list file
    for file in `ls *.$1 2> md5.err`;do
       ((makefile++))
       echo "Making md5 for file $file"
       openssl md5 $file > $file.md5
    done
    rm md5.err
    return 1
}

cleanmd5 () { # clean md5 for *.list file
    for file in `ls *.$1.md5 2> md5.err`;do
       ((cleanfile++))
       echo "Cleaning md5 for file $file"
       rm $file
    done
    rm md5.err
    return 1
}



checkmd5() {  # check md5 of *.list file with *.list.md5
    for file in `ls *.$1 2> md5.err`;do
        echo "Checking  file $file"
        echo "Using md5 file $file.md5"
        ((checkfile++))
        file1=`openssl md5 $file |awk -F " " '{print $2}'`
        #file2=`cut -d " " -f1 $file.md5`
        file3=`grep "$file1" $file.md5 -c`
        
        #echo $file1
        #echo $file2

        #if [ $file1 != $file2 ]
        if [ $file3 != 1 ]
        then
           ((bad++))
           echo "md5 summary BAD"
        else
           ((good++))
           echo "md5 summary GOOD"
        fi
    done
    rm md5.err
    return 1
}

do_recursive_make() {   # recursive execute "makemd5" in this directory
    makemd5 $1
    for filename in `ls`;do
         if [ -d "$filename" ]
         then
             cd $filename
             do_recursive_make $1
             cd ..
         fi
    done
    return 1
}

do_recursive_clean() {   # recursive execute "cleanmd5" in this directory
    cleanmd5 $1
    for filename in `ls`;do
         if [ -d "$filename" ]
         then
             cd $filename
             do_recursive_clean $1
             cd ..
         fi
    done
    return 1
}


do_recursive_check() {   # recursive execute "checkmd5" in this directory

    checkmd5 $1
    for filename in `ls`;do
         if [ -d "$filename" ]
         then
             cd $filename
             do_recursive_check $1
             cd ..
         fi
    done
    return 1
}



if [ $# -ne 2 ]
then
   echo "Usage: ./md5.sh make|clean|check fastq"
   exit
fi

echo "Parameter1: $1"
echo "Parameter2: $2"
if [ $1 == "make" ]
then
   do_recursive_make $2
   echo "Make md5 for $makefile files"
fi

if [ $1 == "clean" ]
then
   do_recursive_clean $2
   echo "Clean md5 for $cleanfile files"
fi


if [ $1 == "check" ]
then
   do_recursive_check $2
   echo "Check md5 for $checkfile files"
   echo "Good: $good"
   echo "Bad:  $bad"
fi

