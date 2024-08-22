#!/bin/perl

#This takes the ancestry blocks and samples every 1MB to give a genotype 
use strict;
use warnings;
use POSIX;

my $marker_spacing = 1000000;
#Load in the minimum range of data per chromosome
my $genome = "/media/drive_5_usb/Childs/texanus_2021/meta/texanus.parentagecovered.bed";
my %chr_start;
my %chr_end;
open GENOME, $genome;
while(<GENOME>){
  chomp;
  my @a = split(/\t/,$_);
  $chr_start{$a[0]} = $a[1];
  $chr_end{$a[0]} = $a[2];
}
close GENOME;

my %state;
my %end;
while(<STDIN>){
  chomp;
  if ($. == 1){next;}
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $start = $a[1];
  my $end = $a[2];
  my $state = $a[3];
  $state{$chr}{$start} = $state;
  $end{$chr}{$start} = $end;
  
}
#try markers 
foreach my $chr (sort keys %chr_start){
  my $marker = 0;
  until($marker > $chr_end{$chr}){
    if ($marker < $chr_start{$chr}){
      #Nothing
    }else{
      foreach my $window_start (sort keys %{$state{$chr}}){
        if (($window_start <= $marker) and ($end{$chr}{$window_start} >= $marker )){
          print "$chr\t$marker\t$state{$chr}{$window_start}\n";
        }
      }
    }
    $marker+=$marker_spacing;
  }
}
