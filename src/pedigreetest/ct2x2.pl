#!/usr/bin/perl -w
use strict;

# 2x2 contingency table
# usage example:
# ct2x2.pl '1:0,0.2' '5: 0.1, 0.5' < datacolumns

my $col_range_1 = shift;	# e.g. '1: 0,0.28' [ col: min,max  whitespace ignored ]
my $col_range_2 = shift;

$col_range_1 =~ s/\s+//g;
$col_range_2 =~ s/\s+//g;

my ($col1, $r1) = split(":", $col_range_1);
my ($col2, $r2) = split(":", $col_range_2);
$col1--;
$col2--;

my ($lb1, $ub1) = split(/[,]/, $r1);
my ($lb2, $ub2) = split(/[,]/, $r2);

my ($in1_in2_count, $in1_out2_count, $out1_in2_count, $out1_out2_count) = (0, 0, 0, 0);
while (<>) {
  my @cols = split(" ", $_);
  my $v1 = $cols[$col1];
  my $v2 = $cols[$col2];
  if ($v1 >= $lb1  and  $v1 < $ub1) { # in 1
    if ($v2 >= $lb2  and  $v2 < $ub2) { # in 2
      $in1_in2_count++;
    } else {			# out 2
      $in1_out2_count++;
    }
  }else{ # out 1
   if ($v2 >= $lb2  and  $v2 < $ub2) { # in 2
      $out1_in2_count++;
    } else {			# out 2
      $out1_out2_count++;
    }
 }
}
my $in1_count = $in1_in2_count + $in1_out2_count;
my $out1_count = $out1_in2_count + $out1_out2_count;
my $in2_count = $in1_in2_count + $out1_in2_count;
my $out2_count = $in1_out2_count + $out1_out2_count;
print "col1: $col1; range 1: $r1 \n";
print "col2: $col2; range 2: $r2 \n\n";

printf("%10s  |%8s %8s  | \n", "", "v2 in", "v2 out");
print("------------------------------------------------\n");
printf("%10s  |%8d %8d  |%8d \n", "v1 in", $in1_in2_count, $in1_out2_count, $in1_count);
printf("%10s  |%8d %8d  |%8d \n", "v1 out", $out1_in2_count, $out1_out2_count, $out1_count);
print("------------------------------------------------\n");
printf("%10s  |%8d %8d  |%8d \n", "v2 totals", $in2_count, $out2_count, $in2_count + $out2_count);

my $N = $in1_in2_count * $out1_out2_count;
my $D = $in1_out2_count * $out1_in2_count;
my $odds_ratio;
if($D > 0){
  $odds_ratio = $N/$D;
  printf("odds ratio: %10.1f \n", $odds_ratio);
}else{
  print "odds ratio: inf \n";
}
