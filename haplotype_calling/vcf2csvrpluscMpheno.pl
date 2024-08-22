#!/bin/perl
use warnings;
use strict;

#This takes a phenotype file, a genetic map file and a vcf file. It outputs a csvr file.
my $min_dp = "3";
my $bc1 = "TRUE"; #Change this to false if it is not a BC1 and if you don't want to make 1/1 genotypes treated as 0/1
my $phenofile = $ARGV[0]; #Tab delimited, first column is sample name
my $mapfile = $ARGV[1]; #Genetic map file
my $snpfile = $ARGV[2]; #Tells whether the alternate or reference allele is from the haplotype of interest. e.g. texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.snpparentage.deb.txt
#Pipe in the vcf file

my %phenotypes;
my %phenohash;
my %samples;
my @samples;
open PHENO, $phenofile;
while(<PHENO>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. == 1){
    foreach my $i (1..$#a){
      $phenotypes{$i} = $a[$i];
    }
  }else{
    my $sample = $a[0];
    $samples{$sample}++;
    push(@samples, $sample);
    foreach my $i (1..$#a){
      $phenohash{$sample}{$phenotypes{$i}} = $a[$i];
    }
  }
}
close PHENO;

my %hash;
open MAP, $mapfile;
while(<MAP>){
  chomp;
  my @a = split(/\t/,$_);
  if ($. ne "1"){
    if ($a[2] ne "NA"){
      my $chrom = $a[0];
      $hash{$chrom}{$a[1]} = $a[2];
    }
  }
}

my %snp_hash;
open SNP, $snpfile;
my %markers;
my %major;
my %minor;
while(<SNP>){
  chomp;
  my @a = split(/\t/,$_);
  my $chr = $a[0];
  my $pos = $a[1];
  my $parent = $a[2];
  my $minor = $a[5];
  my $ref = $a[3];
  my $alt = $a[4];
  $markers{$chr}{$pos} = $minor;
  if ($minor eq  "alt"){
    $major{$chr}{$pos} = $ref;
    $minor{$chr}{$pos} = $alt;
  }else{
    $major{$chr}{$pos} = $alt;
    $minor{$chr}{$pos} = $ref;
  }
}
my %name;
my %reverse_lookup;
my $firstline;
close MAP;
while(<STDIN>){
  chomp;
  my $line = $_;
  my @a = split(/\t/,$line);
  if ($_ =~ m/##/){next;}
  if ($_ =~ m/#/){
    foreach my $i (9..$#a){
      $name{$i} = $a[$i];
      $reverse_lookup{$a[$i]}++;
   }
   foreach my $pheno (sort values %phenotypes){
      unless ($firstline){
        print "$pheno,,";
        $firstline++;
      }else{
        print "\n$pheno,,";
      }
      foreach my $i (9..$#a){
        if (($reverse_lookup{$name{$i}}) and ($samples{$name{$i}})){
          print ",$phenohash{$name{$i}}{$pheno}";
        }
      }
    }
  }else{
    my $chr = $a[0];
    my $chr_numeric = $chr;
    $chr_numeric =~ s/Ha412HOChr//g;
    if ($chr_numeric =~ m/c/){next;} #Skip contigs not in chromosomes.
    my $pos = $a[1];
    my $ref = $a[3];
    my @alt = split(/,/,$a[4]);
    my $bp = $a[1];
    my $previous_site;
    my $before_site;
    my $after_site;
    my $loci_cM;
    foreach my $site (sort  {$a <=> $b} keys %{$hash{$chr}}){
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
    BADSITE:
    next; # Skip site if no cm found.
    FOUNDPOS:
    my $cM_range = $hash{$chr}{$after_site} - $hash{$chr}{$before_site};
    my $bp_range = $after_site - $before_site;
    my $percent_of_range = ($bp - $before_site)/$bp_range;
    $loci_cM = ($percent_of_range * $cM_range) + $hash{$chr}{$before_site};
    
    print "\n${chr}_$pos";
    print ",$chr_numeric,$loci_cM";
   # Section for finding major and minor allele;
    my $major;
    my $minor;
    unless ($markers{$chr}{$pos}){next;} #Only use these markers
    if ($major{$chr}{$pos} eq $ref){
      $major = 0;
    }
    if ($minor{$chr}{$pos} eq $ref){
      $minor = 0;
    }
    foreach my $i (0..$#alt){
      if ($major{$chr}{$pos} eq $alt[$i]){
        $major = $i+1;
      }
      if ($minor{$chr}{$pos} eq $alt[$i]){
        $minor = $i+1;
      }
    }
    foreach my $i (9..$#a){
      unless ($samples{$name{$i}}){next;}
      my @info = split(/:/,$a[$i]);
      my $call = $info[0];
      if ($call eq './.'){ print ",-";next;}
      if ($call eq '.'){print ",-";next;}
      my @genotypes = split(/\//,$call);
      my $minor_count = 0;
      my $major_count = 0;
      foreach my $genotype (@genotypes){
        if ($genotype == $major){
          $major_count++;
        }elsif ($genotype == $minor){
          $minor_count++;
        }
      }
      my $geno = "-";
      if ($minor_count == 2){
        if ($bc1 eq "TRUE"){
          $geno = "H"; #For the BC1 design;
        }else{
          $geno = "B";
        }
      }elsif (($minor_count == 1) and ($major_count == 1)){
          $geno = "H";
      }elsif ($major_count == 2){
          $geno = "A";
      }
      my $dp = $info[1];
      if ($dp eq '.'){$dp = 0;}
      if ($dp < $min_dp){
        $geno = "-";
      }
      print ",$geno";
    }
  }
}

