#!usr/bin/perl
use strict;
use warnings;
use POSIX;

my @keygen = 0;
my @encaps = 0;
my @decaps = 0;
my @keygen_result = 0;
my @encaps_result = 0;
my @decaps_result = 0;

# PARAMS
my $exec = $ARGV[0] . ' > tmp_m';
my $tries = 10000;
my $q1 = 2501;
my $m  = 5001;
my $q3 = 7501; 

for(my $i = 0 ; $i < $tries ; $i = $i +1) {

  if($i % 500 == 0) {
    printf $i . ' ';
  }

  system($exec);
  open(FILE, 'tmp_m');

  while (<FILE>) {
    chomp ;
    if (/Keygen/) {
      @keygen=split(' ');
    }

    if (/Encaps/) {
      @encaps=split(' ');
    }

    if (/Decaps/) {
      @decaps=split(' ');
    }
  }

  $keygen_result[$i] = $keygen[1];
  $encaps_result[$i] = $encaps[1];
  $decaps_result[$i] = $decaps[1];

  close FILE; 
}

@keygen_result = sort { $a <=> $b } @keygen_result;
@encaps_result = sort { $a <=> $b } @encaps_result;
@decaps_result = sort { $a <=> $b } @decaps_result;

print "\n\nKeygen";
printf "\n%d", $keygen_result[$q1];
printf "\n%d", $keygen_result[$m];
printf "\n%d", $keygen_result[$q3];

print "\n\nEncaps";
printf "\n%d", $encaps_result[$q1];
printf "\n%d", $encaps_result[$m];
printf "\n%d", $encaps_result[$q3];

print "\n\nDecaps";
printf "\n%d", $decaps_result[$q1];
printf "\n%d", $decaps_result[$m];
printf "\n%d", $decaps_result[$q3];
printf "\n\n";

system('rm -f tmp_m');

