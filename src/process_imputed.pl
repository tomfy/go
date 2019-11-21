#!/usr/bin/perl -w
use strict;

my $delta = shift // 0.1;


while (<>) {
  if (/^MARKER/) {
    print;
  } else {
    my @cols = split(" ", $_);
    my $seq_id = shift @cols;
    for (my $i = 0; $i < scalar @cols; $i++) {
      my $impval = $cols[$i];
      my $newval = 'X';
      if ($impval < $delta) {
	$newval = 0.0;
      } elsif ($impval > (2.0 - $delta)) {
	$newval = 2.0;
      } elsif ($impval > (1.0 - $delta)  and  $impval < (1.0 + $delta)) {
	$newval = 1.0;
      }
      $cols[$i] = $newval;
    }
    print "$seq_id\t", join("\t", @cols), "\n";
  }
}
