#!/usr/bin/perl

while(glob("*.fas")){
print "/share/raid1/genome/bin/clustalw $_ > log &\n";

}
