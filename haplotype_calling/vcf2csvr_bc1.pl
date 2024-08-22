#!/bin/perl
use warnings;
use strict;
#This takes a vcf file, and converts it to the csvr format of rqtl
#It checks the AF to find the major allele and converts homozygous minor allele to heterozygous because of BC1 design
my %samples;
my $min_dp = 10;
my %name;
my @samples;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/##/){next;}
  if ($_ =~ m/#/){
    print "pheno\t";
    foreach my $i (9..$#a){
        $name{$i} = $a[$i];
        push(@samples,$i);
        print "\t$i";
    }
  }else{
    my $chr = $a[0];
    my $chr_numeric = $chr;
    $chr_numeric =~ s/Ha412HOChr//g;
    if ($chr_numeric =~ m/c/){next;} #Skip contigs not in chromosomes.
    my $pos = $a[1];
    print "\n$chr.$pos\t$chr_numeric";
    my @infos = split(/;/,$a[7]);
    my $major = "ref";
    my $ref_count = 0;
    my $alt_count = 0;
    foreach my $i (@samples){
      my @fields = split(/:/,$a[$i]);
      if ($fields[0] eq "0/0"){
        $ref_count+=2;
      }elsif ($fields[0] eq "0/1"){
        $ref_count++;
        $alt_count++;
      }elsif ($fields[0] eq "1/1"){
        $alt_count+=2;
      }
    }
    my $af = $alt_count/($ref_count + $alt_count);
    if ($af > 0.5){$major = "alt"}
    foreach my $i (@samples){
      my @fields = split(/:/,$a[$i]);
      my $geno = "-";
      if ($fields[0] eq "0/0"){
        if ($major eq "ref"){
          $geno = "A";
        }else{
          $geno = "H";
        }
      }elsif ($fields[0] eq "0/1"){
        $geno = "H";
      }elsif ($fields[0] eq "1/1"){
        if ($major eq "ref"){
       	  $geno = "H"; #For the BC1 design;
        }else{ 
          $geno = "A";
        }
      }
      my $dp = $fields[1];
      if ($dp eq '.'){$dp = 0;}
      if ($dp <= $min_dp){
        $geno = "-";
      }
      print "\t$geno";
    }
  }
}
