#!/usr/bin/perl -w
use strict;

my $file=shift; 


open IN,"$file";
my $first=<IN>;
my $rawreads=(split /,/,$first)[-1];
print $first;
while(<IN>){
    chomp;
    if($_ =~ /%/){
        print $_."\n";
    }else{
        my @tmp=split /,/;
        my $per=($tmp[-1]/$rawreads)*100;
        $per=sprintf("%.2f",$per);
        print "$tmp[0],$per"."%\n";
    }

}close IN;