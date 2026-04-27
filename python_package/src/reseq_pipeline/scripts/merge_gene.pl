#!/usr/bin/perl -w
open IN,$ARGV[0];
$line="genes";
while(<IN>){
	chomp;
	@l=split;
	$line.=",".$l[0];
	$file=$l[0].".gene.xls";
	open I,$file;
	while(<I>){
		next if /^_/;
		@ll=split;
		$s{$ll[0]}.=",".$ll[1];
	}
}
print "$line\n";
for $k (sort keys %s){
	print "$k$s{$k}\n";
}
