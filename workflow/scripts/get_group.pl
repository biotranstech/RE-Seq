#!/usr/bin/perl -w 
#args1	group.xls
#args2	fz.xls
#args3	编辑层面统计结果文件
#args4	编辑层面统计路径

open IN,$ARGV[0];
while(<IN>){
	#chomp;
	@l=split;
	next if /^species/;
	$s{$l[0]}=$l[1];
}

open II,$ARGV[1];
while(<II>){
	chomp;
	@l=split /_vs_/,$_;
	$_=~s/\+/_/;
	print "$_\n";
	$out_file=$ARGV[3]."/$_.csv";
	open O,">$out_file";
	@q=();
	open I,$ARGV[2];
	$headline=readline I;
	$headline=~s/\r//;
	chomp($headline);
	@ll=split /\t/,$headline;
	for ($k=1;$k<=$#ll;$k++){
		if ($s{$ll[$k]} eq $l[0] || $s{$ll[$k]} eq $l[1]){
			push @q,$k;
		}
	}
	$line=$ll[0];
	for $k(@q){$line.=",".$ll[$k];}
	print O "$line\n";
	while(<I>){
		chomp;
		@ll=split /\t/,$_;
		$line=$ll[0];
		for $k(@q){$line.=",".$ll[$k];}
		print O "$line\n";
	}
	close I;
	close O;
}
