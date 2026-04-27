#!/usr/bin/perl -w
open IN,$ARGV[0];
print "species\tSpecies\n";
while(<IN>){
	chomp;
	@l=split;
	print "$l[0]\t$l[3]\n";
}
