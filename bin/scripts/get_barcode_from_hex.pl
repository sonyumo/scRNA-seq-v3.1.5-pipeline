#!/usr/bin/perl -w

use strict;

my $hex=shift;
my $need=shift;
my $out=shift;
open OUT,">$out";
my %hash;

open IN, "$hex";
while(<IN>){
    chomp;
    my @tmp=split /\t/,$_;
    $hash{$tmp[0]}=$tmp[1];

}close IN;


open IN, "$need";
while(<IN>){
    chomp;
    if($hash{$_}){
        print OUT $hash{$_}."\n";
    }
#    else{
#    print "$_\n";
#    }

}close IN;
