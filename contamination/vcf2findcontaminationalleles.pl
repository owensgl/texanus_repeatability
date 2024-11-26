#!/bin/perl
use warnings;
use strict;
#This script takes a vcf file. It is looking for alleles that were not in generation 1 of the texanus cross. It outputs the allele frequency of new alleles for target populations

my %target_location;
my %target_gen;

$target_location{"LBJ"}++;
$target_gen{"7"}++;
$target_gen{"6"}++;
my $min_new_freq = 0.25;

my $list = "/home/owens/working/texanus/texanus.sampleinfo.txt";
my $min_gen1_count = 40;
my $min_genX_count = 40;
open LIST, $list;
my %samples;
my %gen;
my %loc;
while(<LIST>){
  chomp;
  my @a = split(/\t/,$_);
  my $name = $a[0];
  my $gen = $a[2];
  my $loc = $a[5];
  my $type = $a[4];
  if ($type ne "BC"){next;} #Remove wild types
  my $seq = $a[9];
  if (($seq ne "GBS2") and ($seq ne "both")){next;} #Remove GBS1 samples
  $samples{$name}++;
  $gen{$name} = $gen;
  $loc{$name} = $loc;
}


print "chr\tpos\tallele\tfreq";
my %name;
my @samples;
my %total_sites;
my %new_sites;
while(<STDIN>){
  chomp;
  my @a = split(/\t/,$_);
  if ($_ =~ m/##/){next;}
  if ($_ =~ m/#/){
    foreach my $i (9..$#a){
      if ($samples{$a[$i]}){
        $name{$i} = $a[$i];
        push(@samples,$i);
      }
    }
  }else{
    my $gen1_counts;
    my %gen1_alleles;
    my $chr = $a[0];
    my $pos = $a[1];
    my $ref = $a[3];
    my @alts = split(/,/,$a[4]);
    #Check for alleles not present in generation 1
    foreach my $i (@samples){
      if ($a[$i] eq "."){next};
      if ($gen{$name{$i}} ne "1"){next;}
      my @fields = split(/:/,$a[$i]);
      if ($fields[0] eq './.'){next;}
      my @genotypes = split(/\//,$fields[0]);
      foreach my $x (0..1){
	$gen1_alleles{$genotypes[$x]}++;
	$gen1_counts++;
      }
    }
    #Require a minimum of genotypes from generation 1
    unless($gen1_counts){next;}
    if ($gen1_counts < $min_gen1_count){next;}
    my %target_alleles;
    my $target_counts;
    foreach my $i (@samples){
      unless ($target_gen{$gen{$name{$i}}}){next;} #Only look at target generation
      unless ($target_location{$loc{$name{$i}}}){next;} #Only look at target location
      if ($a[$i] eq "."){next};
      my @fields = split(/:/,$a[$i]);
      if ($fields[0] eq './.'){next;}
      my @genotypes = split(/\//,$fields[0]);
      foreach my $x (0..1){
        $target_alleles{$genotypes[$x]}++;
        $target_counts++;
      }
    }
    unless($target_counts){next;}
    unless($target_counts >= $min_genX_count){next;}
    foreach my $allele (sort keys %target_alleles){
      unless ($gen1_alleles{$allele}){
        my $new_freq = $target_alleles{$allele}/$target_counts;
        if ($new_freq < $min_new_freq){next;}
        my $allele_value;
        if ($allele eq 0){
          $allele_value = $ref;
	}else{
	  $allele_value = $alts[($allele-1)];
	}
	print "\n$chr\t$pos\t$allele_value\t$new_freq";
      }
    }
  }
}
