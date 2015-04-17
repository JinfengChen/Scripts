#!/usr/local/bin/perl
 
=head1 NAME
 
    find_bam_insertsize.pl
 
=head1 SYNOPSIS
 
    find_bam_insertsize.pl bam outputdir rpath       
        where bam is the input bam file, 
              outputdir is the output directory for writing output files,
              rpath is the path to the directory where R is installed.
 
=head1 DESCRIPTION
 
    This script takes an input bam file, and finds the mean, median, standard deviation
    and IQR of the insert sizes. The mean, median, standard deviation and IQR are estimated
    by taking the first 1 million read-pairs from the bam file (assuming this won't bias
    the estimates). Of these 1 million read-pairs, only the 'proper' read-pairs are used
    to estimate the mean, median, standard deviation and IQR.
    
=head1 VERSION
  
    Perl script last edited 10-Sept-2012.
 
=head1 CONTACT
 
    alc@sanger.ac.uk (Avril Coghlan)
 
=cut
 
# 
# Perl script find_bam_insertsize.pl
# Written by Avril Coghlan (alc@sanger.ac.uk)
# 10-Sept-12.
# Last edited 25-Jan-2013. 
# SCRIPT SYNOPSIS: find_bam_insertsize: script to estimate the mean, median, standard deviation and IQR of insert sizes in a bam file.
#
#------------------------------------------------------------------#
 
# CHECK IF THERE ARE THE CORRECT NUMBER OF COMMAND-LINE ARGUMENTS:
 
use strict;
use warnings;
 
my $num_args               = $#ARGV + 1;
if ($num_args != 3)
{
    print "Usage of find_bam_insertsize.pl\n\n";
    print "perl find_bam_insertsize.pl <bam> <outputdir> <rpath>\n";
    print "where <bam> is the input bam file,\n";
    print "      <outputdir> is the output directory for writing output files,\n";
    print "      <rpath> is the path to the directory where R is installed\n";
    print "For example, >perl find_bam_insertsize.pl out.sorted.markdup.bam\n";
    print "/lustre/scratch108/parasites/alc/ReaprRepeats/Plasmodium\n";
    print "/software/R-2.15-lenny/bin/\n";
    exit;
}
 
# FIND THE PATH TO THE INPUT BAM FILE:                     
 
my $bam                    = $ARGV[0];
 
# FIND THE DIRECTORY TO USE FOR OUTPUT FILES:      
 
my $outputdir              = $ARGV[1];
 
# FIND THE PATH TO THE DIRECTORY WHERE R IS INSTALLED:
 
my $rpath                  = $ARGV[2];
 
#------------------------------------------------------------------#
 
# TEST SUBROUTINES: 
 
my $PRINT_TEST_DATA        = 0;   # SAYS WHETHER TO PRINT DATA USED DURING TESTING.
&test_print_error;
&test_get_insertsize_stats($outputdir,$rpath);
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
&run_main_program($outputdir,$bam,$rpath);
 
print STDERR "FINISHED.\n";
 
#------------------------------------------------------------------#
 
# RUN THE MAIN PART OF THE CODE:
 
sub run_main_program
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $bam                 = $_[1]; # THE INPUT BAM FILE
   my $rpath               = $_[2]; # PATH TO THE DIRECTORY WHERE R IS INSTALLED
   my $errorcode;                   # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg;                    # RETURNED AS 'none' IF THERE IS NO ERROR. 
   my $output;                      # TEMPORARY OUTPUT FILE FOR STORING THE INSERT SIZES.
   my $median;                      # MEDIAN INSERT SIZE
   my $mean;                        # MEAN INSERT SIZE
   my $IQR;                         # IQR OF INSERT SIZE
   my $sd;                          # STANDARD DEVIATION OF INSERT SIZE
   my $length;                      # SAMPLE SIZE OF INSERT SIZES 
   my $max;                         # MAXIMUM INSERT SIZE 
 
   # READ IN THE INSERT SIZES FROM THE INPUT BAM FILE:
   ($output,$errorcode,$errormsg)  = &read_bam($outputdir,$bam);
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
 
   # GET THE MEAN, MEDIAN, STANDARD DEVIATION AND IQR OF INSERT SIZES:
   ($mean,$median,$sd,$IQR,$length,$max,$errorcode,$errormsg)  = &get_insertsize_stats($output,$outputdir,$rpath);
   print "mean=$mean median=$median sd=$sd IQR=$IQR length=$length max=$max\n";
   if ($errorcode != 0) { ($errorcode,$errormsg) = &print_error($errormsg,$errorcode,0); }
}
 
#------------------------------------------------------------------#
 
# TEST &get_insertsize_stats
 
sub test_get_insertsize_stats
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $rpath               = $_[1]; # PATH TO THE DIRECTORY WHERE R IS INSTALLED
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAMES.
   my $file;                        # TEMPORARY FILE TO PUT NUMBERS INTO.
   my $median;                      # MEDIAN INSERT SIZE
   my $mean;                        # MEAN INSERT SIZE
   my $IQR;                         # IQR OF INSERT SIZE
   my $sd;                          # STANDARD DEVIATION OF INSERT SIZE
   my $length;                      # SAMPLE SIZE OF INSERT SIZES 
   my $max;                         # MAXIMUM INSERT SIZE
   my $errorcode;                   # RETURNED AS 0 BY A FUNCTION IF THERE WAS NO ERROR
   my $errormsg;                    # RETURNED AS 'none' BY A FUNCTION IF THERE WAS NO ERROR
 
   $random_number            = rand();
   $file                     = $outputdir."/tmp".$random_number;
   open(FILE,">$file") || die "ERROR: test_get_insertsize_stats: cannot open $file\n";
   print FILE "1232\n";
   print FILE "1263\n";
   print FILE "2321\n";
   print FILE "523\n";
   print FILE "233\n";
   print FILE "432\n";
   close(FILE); 
   ($median,$mean,$IQR,$sd,$length,$max,$errorcode,$errormsg)    = &get_insertsize_stats($file,$outputdir,$rpath);
   if ($errorcode != 0 || $median != 877.5 || $mean != 1000.667 || $IQR != 800.5 || $sd != 775.4319 || $length != 6)
   {
      print STDERR "ERROR: test_get_insertsize_stats: failed test1\n";
      exit;
   }
}
 
#------------------------------------------------------------------#
 
# GET THE MEAN, MEDIAN, STANDARD DEVIATION AND IQR OF INSERT SIZES:
 
sub get_insertsize_stats
{
   my $output              = $_[0]; # TEMPORARY FILE WITH THE INSERT SIZES
   my $outputdir           = $_[1]; # DIRECTORY FOR WRITING OUTPUT INTO 
   my $rpath               = $_[2]; # PATH TO THE DIRECTORY WHERE R IS INSTALLED
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR
   my $errormsg            = 'none';# RETURNED AS 'none' IF THERE IS NO ERROR
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
   my $script;                      # TEMPORARY FILE FOR R SCRIPT 
   my $median;                      # MEDIAN INSERT SIZE
   my $mean;                        # MEAN INSERT SIZE
   my $IQR;                         # IQR OF INSERT SIZE
   my $sd;                          # STANDARD DEVIATION OF INSERT SIZE
   my $length;                      # SAMPLE SIZE OF INSERT SIZES 
   my $max;                         # MAXIMUM INSERT SIZE
   my $output2;                     # OUTPUT FILE TO WRITE TO IN R 
   my $line;                        # 
   my @temp;                        # 
 
   $random_number            = rand();
   $script                   = $outputdir."/tmp".$random_number;
   $random_number            = rand();
   $output2                  = $outputdir."/tmp".$random_number;
   open(SCRIPT,">$script") || die "ERROR: get_insertsize_stats: cannot open $script\n";
   print SCRIPT "#!$rpath/Rscript\n";
   print SCRIPT "MyData <- read.table(\"$output\",header=FALSE)\n";
   print SCRIPT "insert <- MyData\$V1\n";
   print SCRIPT "mymedian <- median(insert)\n";
   print SCRIPT "mymean <- mean(insert)\n";
   print SCRIPT "sd <- sd(insert)\n";
   print SCRIPT "iqr <- IQR(insert)\n";
   print SCRIPT "max <- max(insert)\n";
   print SCRIPT "len <- length(insert)\n";
   print SCRIPT "mydata <- c(mymedian, mymean, sd, iqr, max, len)\n";
   print SCRIPT "write(mydata,file=\"$output2\",ncolumns=6)\n";
   close(SCRIPT);
   system "chmod +x $script";
   system "$script";
   system "rm -f $script";
   # READ IN THE MEDIAN, MEAN, IQR, STANDARD DEVIATION, SAMPLE SIZE AND MAXIMUM INSERT SIZE:
   open(OUTPUT2,"$output2") || die "ERROR: get_insertsize_stats: cannot open $output2\n";
   while(<OUTPUT2>)
   {
      $line                = $_;
      @temp                = split(/\s+/,$line);
      $median              = $temp[0];
      $mean                = $temp[1];
      $sd                  = $temp[2];
      $IQR                 = $temp[3];
      $max                 = $temp[4];
      $length              = $temp[5];
   }
   close(OUTPUT2);
   # REMOVE THE TEMPORARY DATA FILE:
   system "rm -f $output"; 
   system "rm -f $output2";
 
   return($median,$mean,$IQR,$sd,$length,$max,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# READ IN THE INSERT SIZES FROM THE INPUT BAM FILE:
 
sub read_bam            
{
   my $outputdir           = $_[0]; # DIRECTORY TO PUT OUTPUT FILES IN.
   my $bam                 = $_[1]; # THE INPUT BAM FILE
   my $errorcode           = 0;     # RETURNED AS 0 IF THERE IS NO ERROR.
   my $errormsg            = "none";# RETURNED AS 'none' IF THERE IS NO ERROR.
   my $output;                      # TEMPORARY OUTPUT FILE FOR PUTTING THE INSERT SIZES INTO 
   my $random_number;               # RANDOM NUMBER TO USE IN TEMPORARY FILE NAME
   my $output2;                     # TEMPORARY OUTPUT FILE FOR PUTTING THE INSERT SIZES INTO
   my $line;                        # 
   my @temp;                        # 
   my $insertsize;                  # INSERT SIZE
   my $flag;                        # THE SAMTOOLS BIT FLAG 
 
   # OPEN A TEMPORARY FILE FOR PUTTING THE INSERT SIZES INTO:   
   $random_number            = rand();
   $output                   = $outputdir."/tmp".$random_number;
   open(OUTPUT,">$output") || die "ERROR: read_bam: cannot open $output\n";
   
   # GET THE INSERT SIZES FROM THE BAM FILE AND PUT THEM INTO THE TEMPORARY OUTPUT FILE:
   system "samtools view $bam | cut -f2,9 | head -1000000 > $output"; # USE A SAMPLE OF 1000,000 READ-PAIRS
   close(OUTPUT);
 
   # NOW JUST EXTRACT THE INSERT SIZES WITH GOOD BIT FLAGS
   # (FROM http://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/
   # THESE ARE 99, 147, 83 OR 163):
   $random_number            = rand();
   $output2                  = $outputdir."/tmp".$random_number;
   open(OUTPUT2,">$output2") || die "ERROR: read_bam: cannot open $output2\n";
   open(OUTPUT,"$output") || die "ERROR: read_bam: cannot open $output\n";
   while(<OUTPUT>)
   { 
      $line                  = $_;
      chomp $line;
      @temp                  = split(/\s+/,$line);
      $flag                  = $temp[0];
      $insertsize            = $temp[1];
      # ONLY 99, 147, 83 OR 163 ARE PROPERLY MAPPED READ PAIRS WITHIN A DEFINED INSERT SIZE
      if ($flag == 99 || $flag == 147 || $flag == 83 || $flag == 163)
      {
         if ($flag == 99 || $flag == 163)
         {
            if ($insertsize < 0) 
            {
               $errormsg     = "ERROR: read_bam: insertsize $insertsize but flag $flag\n";
               $errorcode    = 1; # ERRORCODE=1
               return($output2,$errorcode,$errormsg);
            }
            print OUTPUT2 "$insertsize\n";
         }
         elsif ($flag == 147 || $flag == 83)
         {
            if ($insertsize > 0)
            {
               $errormsg     = "ERROR: read_bam: insertsize $insertsize but flag $flag\n";
               $errorcode    = 2; # ERRORCODE=2
               return($output2,$errorcode,$errormsg);
            } 
            $insertsize      = -$insertsize;
            print OUTPUT2 "$insertsize\n";
         }
         else
         {
            $errormsg        = "ERROR: read_bam: flag $flag\n";
            $errorcode       = 3; # ERRORCODE=3
            return($output2,$errorcode,$errormsg);
         }
      }
   }
   close(OUTPUT); 
   close(OUTPUT2);
 
   return($output2,$errorcode,$errormsg);
}
 
#------------------------------------------------------------------#
 
# TEST &print_error
 
sub test_print_error
{
   my $errormsg;                    # RETURNED AS 'none' FROM A FUNCTION IF THERE WAS NO ERROR
   my $errorcode;                   # RETURNED AS 0 FROM A FUNCTION IF THERE WAS NO ERROR
 
   ($errormsg,$errorcode)  = &print_error(45,45,1);
   if ($errorcode != 12) { print STDERR "ERROR: test_print_error: failed test1\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('My error message','My error message',1);
   if ($errorcode != 11) { print STDERR "ERROR: test_print_error: failed test2\n"; exit;}
 
   ($errormsg,$errorcode)  = &print_error('none',45,1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test3\n"; exit;} 
 
   ($errormsg,$errorcode)  = &print_error('My error message', 0, 1);
   if ($errorcode != 13) { print STDERR "ERROR: test_print_error: failed test4\n"; exit;}
}
 
#------------------------------------------------------------------#
 
# PRINT OUT AN ERROR MESSAGE AND EXIT.
 
sub print_error
{
   my $errormsg            = $_[0]; # THIS SHOULD BE NOT 'none' IF AN ERROR OCCURRED.
   my $errorcode           = $_[1]; # THIS SHOULD NOT BE 0 IF AN ERROR OCCURRED.
   my $called_from_test    = $_[2]; # SAYS WHETHER THIS WAS CALLED FROM test_print_error OR NOT
 
   if ($errorcode =~ /[A-Z]/ || $errorcode =~ /[a-z]/) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 11; $errormsg = "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; # ERRORCODE=11
         return($errormsg,$errorcode);
      }
      else 
      { 
         print STDERR "ERROR: print_error: the errorcode is $errorcode, should be a number.\n"; 
         exit;
      }
   }
 
   if (!($errormsg =~ /[A-Z]/ || $errormsg =~ /[a-z]/)) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 12; $errormsg = "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; # ERRORCODE=12
         return($errormsg,$errorcode);
      }
      else
      {
         print STDERR "ERROR: print_error: the errormessage $errormsg does not seem to contain text.\n"; 
         exit;
      }
   }
 
   if    ($errormsg eq 'none' || $errorcode == 0) 
   { 
      if ($called_from_test == 1)
      {
         $errorcode = 13; $errormsg = "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; # ERRORCODE=13
         return($errormsg,$errorcode);
      }
      else 
      {
         print STDERR "ERROR: print_error: errormsg $errormsg, errorcode $errorcode.\n"; 
         exit;
      }
   }
   else                                           
   { 
      print STDERR "$errormsg"; 
      exit;                                                      
   } 
 
   return($errormsg,$errorcode);
}
 
#------------------------------------------------------------------#


