#!/usr/bin/perl -w
use strict;

my $gfile = shift; # genotypes file (matrix format, rows: accessions, cols: markers)
my $pfile = shift; # pedigree file,

open my $fhg, "<", $gfile or die "Couldn't open $gfile for reading.\n";
my $firstline = <$fhg>;
my $n_markers = (scalar split(" ", $firstline) ) - 1;
my %id_genotypes;
while (my $line = <$fhg>) {
  my @fields = split(" ", $line);
  my $id = shift @fields;
  $id_genotypes{$id} = join('', @fields);
}
close $fhg;

open my $fhp, "<", "$pfile" or die "Couldn't open $pfile for reading.\n";
my $pfirstline = <$fhp>;
while (my $line = <$fhp>) {
  my @fields = split(" ", $line);
  my ($acc_id, $mat_id, $pat_id) = @fields[-3, -2, -1]; # take the last 3 fields
  if($acc_id eq 'NA'  or  $mat_id eq 'NA'  or   $pat_id eq 'NA'){ # missing information - skip.
    warn "skipping. $acc_id  $mat_id  $pat_id \n";
  }else{
  my $a_gts = $id_genotypes{$acc_id} // undef;
  my $m_gts = $id_genotypes{$mat_id} // undef;
  my $p_gts = $id_genotypes{$pat_id} // undef;
  if (defined $a_gts  and  defined $m_gts  and  defined $p_gts) {
    my ($goods, $bads) = good_bad_counts($a_gts, $m_gts, $p_gts);
    print "$acc_id  $mat_id  $pat_id    $goods  $bads  ",  $bads/($goods + $bads), "\n";
  }else{
    warn "One of the 3 ids (accession of parents) not present in genotypes matrix file.\n"
  }
}


}

sub good_bad_counts{
  my $acc_gts = shift;		# string
  my $mat_gts = shift;
  my $pat_gts = shift;

  my ($good_count, $bad_count) = (0, 0);

  for (my $i = 0; $i < length $acc_gts; $i++) {
    my $a_gt = substr($acc_gts, $i, 1);
    my $m_gt = substr($mat_gts, $i, 1);
    my $p_gt = substr($pat_gts, $i, 1);
    if ($m_gt == 1  and  $p_gt == 1) {
      # $a_gt can be 0, 1, or 2. So doesn't add any info. Skip for now.
    } elsif($m_gt == $p_gt) { # i.e. either 0,0 or 2,2
      if (is_possible($a_gt, $m_gt, $p_gt) == 1) {
	$good_count++;
      } else {
	$bad_count++;
      }
    }
  }
  return ($good_count, $bad_count);
}

# returns 1 if $a_gt could result from cross of $m_gt and $p_gt, 0 otherwise.
# all 3 arguments assumed to be 0, 1, or 2.
sub is_possible{
  my $a_gt = shift;
  my $m_gt = shift;
  my $p_gt = shift;

  my ($g_0022, $b_0022, $g_0110, $b_0110, $g_1221, $b_1221, $g_0220, $b_0220, $g_11) = (0, 0, 0, 0, 0, 0, 0, 0, 0);

  if ($m_gt == 0) {
    if ($p_gt == 0) {
      return ($a_gt == 0)? 1 : 0;
    } elsif ($p_gt == 1) {
      return ($a_gt == 2)? 0 : 1;
    } else {			# $p_gt == 2
      return ($a_gt == 1)? 1 : 0;
    }
  } elsif ($m_gt == 1) {
    if ($p_gt == 0) {
      return ($a_gt == 2)? 0 : 1;
    } elsif ($p_gt == 1) {
      return 1;
    } else {			# $p_gt == 2
      return ($a_gt == 0)? 0 : 1;
    }
      
  } else {			# $m_gt == 2
    if ($p_gt == 0) {
      return ($a_gt == 1)? 1 : 0;
    } elsif ($p_gt == 1) {
      return ($a_gt == 0)? 0 : 1;
    } else {			# $p_gt == 2
      return ($a_gt == 2)? 1 : 0;
    }
  }
  die "Shouldn't get here.\n";
}
      
