#!/usr/bin/perl -w
use strict;

# for each marker, get the number of 0s, 1s, 2s, Xs.
# read matrix file from stdin.

my $first_line = <>;
my @markers = split(" ", $first_line);

my %marker_count = ();
my $x = shift @markers;
for (@markers) {
  $marker_count{$_} = [0, 0, 0, 0];
}
while (my $line = <>) {
  my @gts = split(" ", $line);
  my $id = shift @gts;
  while (my ($i, $gt) = each @gts) {
    my $m = $markers[$i];
    if ($gt eq 'X') {
      $marker_count{$m}->[3]++;
    } else {
      $marker_count{$m}->[$gt]++;
    }
  }
}

while (my ($m, $gtc) = each %marker_count) {
  my $n012 = $gtc->[0] + $gtc->[1] + $gtc->[2];
  printf("%16s  %6d %6d %6d %6d %7.4f %7.4f %7.4f \n", $m, $gtc->[0], $gtc->[1], $gtc->[2], $gtc->[3],
	 $gtc->[0]/$n012, $gtc->[1]/$n012, $gtc->[2]/$n012);
}

