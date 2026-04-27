#!/usr/bin/perl

use strict;
use warnings;

my $input_file = $ARGV[0] . '.hardfiltered.biallelic.vcf';
my $output_file = $ARGV[0] . '.Reads.csv';
## 打开 VCF 文件进行读取
open(my $vcf_fh, '<', $input_file) or die "Cannot open input file: $!";
## 打开输出 CSV 文件进行写入
open(my $csv_fh, '>', $output_file) or die "Cannot create output file: $!";

## 打印 CSV 文件的表头
print $csv_fh "Sample,Chr,Pos,Ref,Alt,DP,ReadsCount,RefReadsCount,MutationReadsCount,MutationFrequency\n";

## 逐行读取 VCF 文件
while (my $line = <$vcf_fh>) {
    chomp $line;
    next if $line =~ /^#/;  # 跳过注释行
     # 解析 VCF 行
    my @fields = split /\t/, $line;
    my $chrom = $fields[0];
    my $pos = $fields[1];
    my $ref = $fields[3];
    my $alt = $fields[4];
    my $info = $fields[7];
    my $readeds = $fields[9];
    # 提取总 Reads 数量
    my ($DP) = $info =~ /DP=(\d+)/;

    # 提取突变和参考 Reads 数量
    my ($ref_reads, $mutation_reads, $reads_count) = (0, 0, 0);
    if ($readeds =~ /:(\d+),(\d+):(\d+)/) {
	($ref_reads, $mutation_reads, $reads_count) = ($1, $2, $3);
    } elsif ($readeds =~ /(\d+),(\d+)/) {
	($ref_reads, $mutation_reads) = ($1, $2);
	$reads_count = $ref_reads + $mutation_reads;
    }

    # 计算突变频率
    my $mutation_frequency;
    if (defined $mutation_reads && $reads_count > 0) {
        $mutation_frequency = $mutation_reads / $reads_count;
    } else {
        $mutation_frequency = 0;
    }
    #print "Reads: $reads\n";
    # 将数据写入 CSV 文件
    print $csv_fh "$ARGV[0],$chrom,$pos,$ref,$alt,$DP,$reads_count,$ref_reads,$mutation_reads,$mutation_frequency\n";
}
# 关闭文件句柄
close($vcf_fh);
close($csv_fh);
print "Output CSV file generated successfully.\n";
                                                                        
