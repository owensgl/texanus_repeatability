#!/bin/perl
#This script takes a csvr and outputs the allele frequency and count for each site
use strict;
use warnings;
print "marker\tpercent_A\tgenotyped";
while(<STDIN>){
  chomp;
  unless ($_ =~ /^Ha412/){
    next;
  }
  my @a = split(/,/,$_);
  my $marker = $a[0];
  my $n_count = 0;
  my $A_count = 0;
  foreach my $i (3..$#a){
    if ($a[$i] eq "-"){next;}
    $n_count++;
    if ($a[$i] eq "A"){
      $A_count++;
    }
  }
  my $percent = "NA";
  if ($n_count > 0){
    $percent = $A_count/$n_count;
  }
  print "\n$marker\t$percent\t$n_count";
}

