#!/bin/perl
#This script takes a set of ~25% snps from the BC1 population, converts false indels to SNPs and then compares them to the debilis vcf to find if the allele is common in debilis. It skips real indels;
use strict;
use warnings;

my $bc1_vcf = $ARGV[0];

open(BC1, "cat $bc1_vcf |");

print "chr\tpos\toriginal_pos\tref\talt";
print "\tminor\taf";
my %sites;
my %alleles;
while(<BC1>){
	if ($_ =~ /^#/){next;}
	my @a = split(/\t/,$_);
	my $chr = $a[0];
	my $pos = $a[1];
        my $original_pos = $pos;
	my $ref = $a[3];
	my $alt = $a[4];
	my $indel;
	my $ref_snp;
	my $alt_snp;
	my $snp_site;
	#Check to see which allele is more common
	my $ref_count = 0;
	my $alt_count = 0;
	foreach my $i (9..$#a){
		my @fields = split(/:/,$a[$i]);
		if ($fields[0] eq "0/0"){
			$ref_count+=2;
		}elsif($fields[0] eq "0/1"){
			$ref_count+=1;
			$alt_count+=1;
		}elsif($fields[0] eq "1/1"){
			$alt_count+=2;
		}
	}
	my $af = $alt_count / ($ref_count + $alt_count);
	my $minor = "alt";
	if ($af > 0.5){
		$minor = "ref";
	}
	unless($ref_snp){
		$ref_snp = $ref;
	}
	unless($alt_snp){
		$alt_snp = $alt;
	}
	unless($snp_site){
		$snp_site = $pos;
	}
	print "\n$chr\t$snp_site\t$original_pos\t$ref_snp\t$alt_snp";
	print "\t$minor\t$af";
	
}
