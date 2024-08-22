#!/bin/perl

use strict;
use warnings;
#This takes a vcf and makes a new file with cm locations for each position in the vcf.
#Pipe in a vcf
my $map = "/media/drive_5_usb/Fuchs/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt";

my %cm_hash;
open MAP, $map;
while(<MAP>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. != 1){
    my $chrom = $a[0];
    $cm_hash{$chrom}{$a[1]} = $a[2];
  }
}

close MAP;

print "chr\tpos\tcm";
while(<STDIN>){
  chomp;
  my %geno;
  my @a = split(/\t/,$_);
  if ($_ =~ m/^##/){next;}
  if ($_ =~ m/^#/){
  }else{
    my $chrom = $a[0];
    my $bp = $a[1];
    my $previous_site;
    my $before_site;
    my $after_site;
    my $loci_cM;
    foreach my $site (sort  {$a <=> $b} keys %{$cm_hash{$chrom}}){
      if ($site > $bp){
        if ($previous_site){
          $before_site = $previous_site;
          $after_site = $site;
          goto FOUNDPOS;
        }else{
          $loci_cM = "NA";
          goto BADSITE;
        }
      }
      $previous_site = $site;
    }
    $loci_cM = "NA";
    goto BADSITE;
    FOUNDPOS:
    my $cM_range = $cm_hash{$chrom}{$after_site} - $cm_hash{$chrom}{$before_site};
    my $bp_range = $after_site - $before_site;
    my $percent_of_range = ($bp - $before_site)/$bp_range;
    $loci_cM = ($percent_of_range * $cM_range) + $cm_hash{$chrom}{$before_site};
    
    BADSITE:
    print "\n";
    foreach my $i(0..1){
      print "$a[$i]\t";
    }
    print "$loci_cM";
  }
}
