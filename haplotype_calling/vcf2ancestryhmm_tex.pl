#!/bin/perl

use strict;
use warnings;
#This takes a vcf and converts it to a ancestryhmm file. We're using 95% confidence for parent specific markers to account for outside gene flow adding alleles. It uses a cm map from annuus. 
#subset by generation and population before running.  
my $allele_info = $ARGV[0]; #texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.snpparentage.deb.txt
#pipe in the vcf.
my $map = $ARGV[1]; #texanus.freebayes.5dp.parentagesnps.cmlocations.txt #Created using vcf2cmlocations.pl
my $output_prefix = $ARGV[2];
my $sample_order = $ARGV[3]; #Just a list of samples, one column in the order you want them output

my @sample_order;
open SAMPLE, $sample_order;
while(<SAMPLE>){
  chomp;
  push(@sample_order, $_);
}
close SAMPLE;

my %cm_hash;
open MAP, $map;
while(<MAP>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $cm = $a[2];
  $cm_hash{$chr}{$pos} = $cm;
}
close MAP;

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

my %data;
my $target_counts;
my $total_counts;
my %sample;
my $previous_cm = 0;
my %cm_gap_hash;
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
    my $cm = $cm_hash{$chr}{$pos};
    $cm_gap_hash{$chr}{$pos} = $cm - $previous_cm;
    $previous_cm = $cm;
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
        $data{$chr}{$pos}{$sample{$i}} = "0\t0";
        next;
      }else{
        my @infos = split(/:/,$a[$i]);
        if (($infos[0] eq '.') or ($infos[0] eq './.')){
          $data{$chr}{$pos}{$sample{$i}} = "0\t0";
          next;
        }
        my @genos = split(/\//,$infos[0]);
        my $total_dp = $infos[1];
        my @specific_dp = split(/,/,$infos[2]);
        my $target_dp = $specific_dp[$target_allele];
        my $off_target_dp = $total_dp - $target_dp;
        $data{$chr}{$pos}{$sample{$i}} = "$target_dp\t$off_target_dp";
        foreach my $j (0..1){
          if ($genos[$j] eq $target_allele){
            $target_counts++;
            $total_counts++;
          }else{
            $total_counts++;
          }
        }
      }  
    }
  }
}

open (OUT1, '>',"$output_prefix.ancestry_proportion.txt");
open (OUT2, '>',"$output_prefix.ancestry_hmminput.txt");
my $proportions = $target_counts/$total_counts;
print OUT1 "Ancestry_proportion $proportions";

foreach my $chr (sort keys %data){
  foreach my $pos (sort  {$a <=> $b} keys %{$data{$chr}} ){
    print OUT2 "$chr\t$pos\t100\t0\t5\t95\t$cm_gap_hash{$chr}{$pos}";
    foreach my $sample (@sample_order){
     print OUT2 "\t$data{$chr}{$pos}{$sample}";  
    }
    print OUT2 "\n";
  }
}


