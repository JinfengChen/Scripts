#!/usr/bin/perl

print &get_time()."\n";

sub get_time {
   ($sec,$min,$hour,$day,$mon,$year,$weekday,$yeardate,$savinglightday)
     = (localtime(time));
 
   $sec  = ($sec < 10)? "0$sec":$sec;
   $min  = ($min < 10)? "0$min":$min;
   $hour = ($hour < 10)? "0$hour":$hour;
   $day  = ($day < 10)? "0$day":$day;
   $mon  = ($mon < 9)? "0".($mon+1):($mon+1);
   $year += 1900;
   
   return "$year-$mon-$day $hour:$min:$sec.00";
   
}

