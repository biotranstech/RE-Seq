#!/usr/bin/perl -w
open IN,$ARGV[0];
while(<IN>){
	chomp;
	@l=split;
	if (/^#/){print "$_\n";}
	else{
		if(/^#/ || $l[6] eq "PASS" || $l[6] eq "MQRankSum-12.5;ReadPosRankSum-8" ||$l[6] eq "MQRankSum-12.5" ||$l[6] eq "ReadPosRankSum-8"){
		if ($l[4] !~ /[a-zA-Z]*,[a-zA-Z]*/){print "$_\n";}
		}
	}
}
