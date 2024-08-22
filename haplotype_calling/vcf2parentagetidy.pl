#!/bin/perl

use strict;
use warnings;
#This takes a vcf and converts it to a ancestryhmm file. We're using 95% confidence for parent specific markers to account for outside gene flow adding alleles. It uses a cm map from annuus.
#subset by generation and population before running.
my $allele_info = $ARGV[0]; #texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.snpparentage.deb.txt
#pipe in the vcf.


my %allele_info;
open ALLELE, $allele_info;
while(<ALLELE>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $ref = $a[3];
  my $alt = $a[4];
  my $type = $a[5];
  my $target_allele;
  if ($type eq "alt"){
    $target_allele = $alt;
  }else{
    $target_allele = $ref;
  }
  $allele_info{$chr}{$pos} = $target_allele;
}
close ALLELE;
my %pop;

my %sample;
print "chr\tpos\tsample\tgenotype_count\n";
while(<STDIN>){
  chomp;
  my %geno;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
    foreach my $i (9..$#a){
      $sample{$i} = $a[$i];
    }
  }else{
    my $chr = $a[0];
    my $pos = $a[1];
    unless($allele_info{$chr}{$pos}){next;}
    my $ref = $a[3];
    my @alts = split(/,/,$a[4]);
    my $target_allele;
    if ($ref eq $allele_info{$chr}{$pos}){
      $target_allele = 0;
    }else{
      foreach my $i (0..$#alts){
        if ($alts[$i] eq $allele_info{$chr}{$pos}){
          $target_allele = $i+1;
        }
      }
    }

    foreach my $i (9..$#a){
      if (($a[$i] eq '.') or ($a[$i] eq './.')){
        print "$chr\t$pos\t$sample{$i}\tNA\n";
      }else{
        my @infos = split(/:/,$a[$i]);
        if (($infos[0] eq '.') or ($infos[0] eq './.')){
          print "$chr\t$pos\t$sample{$i}\tNA\n";
          next;
        }
        my @genos = split(/\//,$infos[0]);
        my $target_counts =0;
        foreach my $j (0..1){
          if ($genos[$j] eq $target_allele){
            $target_counts++;
          }else{
          }
        }
        print "$chr\t$pos\t$sample{$i}\t$target_counts\n";
      }
    }
  }
}

