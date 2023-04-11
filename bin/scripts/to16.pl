#!/usr/bin/perl -w
use strict;
my $input=shift;



open IN,"$input";
while(my $in=<IN>){
    chomp$in;
    my @tmp=split /\t/,$in;
#    print "$tmp[0]\t$tmp[1]\t";
    $tmp[0] =~ tr/A/0/;
    $tmp[0] =~ tr/C/1/;
    $tmp[0] =~ tr/G/2/;
    $tmp[0] =~ tr/T/3/;
#    print "$tmp[0]\t";
    my @atcg=split //,$tmp[0];
    my $bit10=$atcg[19]*(4**0) + $atcg[18]*(4**1) + $atcg[17]*(4**2)+ $atcg[16]*(4**3)+ $atcg[15]*(4**4)+ $atcg[14]*(4**5)+ $atcg[13]*(4**6)+ $atcg[12]*(4**7)+ $atcg[11]*(4**8)+ $atcg[10]*(4**9)+ $atcg[9]*(4**10)+ $atcg[8]*(4**11)+ $atcg[7]*(4**12)+ $atcg[6]*(4**13)+ $atcg[5]*(4**14)+ $atcg[4]*(4**15)+ $atcg[3]*(4**16)+ $atcg[2]*(4**17)+ $atcg[1]*(4**18)+ $atcg[0]*(4**19);
    my $six=sprintf("%X",$bit10);
    if(length($six) == 9){
        print "0".$six."\t$tmp[1]\n";
    }else{
        print "$six\t$tmp[1]\n";
    }

}close IN;
