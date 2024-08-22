#!/bin/perl

#This takes the posterior probability from ancestry HMM and turns it into blocks.
#Takes boundaries as the midpoint between switches.
use strict;
use warnings;
use POSIX;

print "chr\tstart\tend\tstate";
my $window_chr = "NA";
my $window_start;
my $last_position;
my $window_state;
while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $current_state;
  if ($a[3] > 0.5){
    $current_state = "01";
  }elsif ($a[2] > 0.5){
    $current_state = "00";
  }elsif ($a[4] > 0.5){
    $current_state = "11";
  }else{
    $current_state = "NA";
  }
  #It's a new chromosome
  if ($window_chr ne $chr){
    #print previous window
    if ($window_start){
      print "\n$window_chr\t$window_start\t$last_position\t$window_state";
    }
    #Start new window
    $window_chr = $chr;
    $window_state = $current_state;
    $window_start = $pos;
    
  }elsif ($current_state ne $window_state){
    #If its switched window states
    my $gap_position = floor(($last_position + $pos) / 2);
    print "\n$window_chr\t$window_start\t$gap_position\t$window_state";
    $window_start = $gap_position+1;
    $window_state = $current_state;
  }
  $last_position = $pos;
}

print "\n$window_chr\t$window_start\t$last_position\t$window_state";
