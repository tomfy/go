#!/usr/bin/perl -w
use strict;

my $gfile = shift; # genotypes file (matrix format, rows: accessions, cols: markers)
my $pfile = shift; # pedigree file,
my $p = shift // 0.9; # 90th %ile by default

open my $fhg, "<", $gfile or die "Couldn't open $gfile for reading.\n";
my %id_genotypes;
my $firstline = <$fhg>;
my @marker_ids = split(" ", $firstline);
my $mid = shift @marker_ids;
die "First line of genotype file should start with MARKER. Instead: $mid \n" if ($mid ne "MARKER");
my $n_markers = scalar @marker_ids;

while (my $line = <$fhg>) {
  my @fields = split(" ", $line);
  my $acc_id = shift @fields;
  $id_genotypes{$acc_id} = join('', @fields);
}
close $fhg;
my @ids = keys %id_genotypes;
my $n_ids = scalar @ids;

my @n00s = (0) x $n_markers;
my @n001s = (0) x $n_markers;
my @n002s = (0) x $n_markers;

my @n01s = (0) x $n_markers;
my @n012s = (0) x $n_markers;

my @n02s = (0) x $n_markers;
my @n020s = (0) x $n_markers;
my @n022s = (0) x $n_markers;

my @n12s = (0) x $n_markers;
my @n120s = (0) x $n_markers;

my @n22s = (0) x $n_markers;
my @n220s = (0) x $n_markers;
my @n221s = (0) x $n_markers;

my $nX;

open my $fhp, "<", "$pfile" or die "Couldn't open $pfile for reading.\n";
my $pfirstline = <$fhp>;
while (my $line = <$fhp>) {
  my @fields = split(" ", $line);
  my ($acc_id, $mat_id, $pat_id) = @fields[-3, -2, -1]; # take the last 3 fields
  if ($acc_id eq 'NA'  or  $mat_id eq 'NA'  or   $pat_id eq 'NA') { # missing information - skip.
    warn "skipping. $acc_id  $mat_id  $pat_id \n";
  } else {
    my $a_gts = $id_genotypes{$acc_id} // undef;
    my $m_gts = $id_genotypes{$mat_id} // undef;
    my $p_gts = $id_genotypes{$pat_id} // undef;
    # my $rid1 = $ids[rand($n_ids)];
    # my $r1_gts = $id_genotypes{$rid1};
    # my $r2_gts = $id_genotypes{$ids[rand($n_ids)]};
    if (defined $a_gts  and  defined $m_gts  and  defined $p_gts) {
      # score rel. to parents specified in the pedigree.
      #   my ($goods, $bads) = good_bad_counts($a_gts, $m_gts, $p_gts);
      #   my $score = $bads/($goods + $bads);
      #    print STDERR "Score: $score \n"; 

      my $nX = nxyzs($a_gts, $m_gts, $p_gts,
		     \@n00s, \@n001s, \@n002s,
		     \@n01s, \@n012s,
		     \@n02s, \@n020s, \@n022s,
		     \@n12s, \@n120s,
		     \@n22s, \@n220s, \@n221s);

    } else {
      warn "One of the 3 ids (accession or parents) not present in genotypes matrix file.\n"
    }
  }

}
my @x001221s = ();
my @x002220s = ();
my @x012120s = ();
my @x020022s = ();

while (my ($i, $mid) = each @marker_ids) {

  my $n00 = $n00s[$i];
  my $x001 = ($n00 > 0)? $n001s[$i]/$n00 : '-';
  my $x002 = ($n00 > 0)? $n002s[$i]/$n00 : '-';

  my $n01 = $n01s[$i];
  my $x012 = ($n01 > 0)? $n012s[$i]/$n01 : '-';

  my $n02 = $n02s[$i];
  my $x020 = ($n02 > 0)? $n020s[$i]/$n02 : '-';
  my $x022 = ($n02 > 0)? $n022s[$i]/$n02 : '-';
  
  my $n12 = $n12s[$i];
  my $x120 = ($n12 > 0)? $n120s[$i]/$n12 : '-';

  my $n22 = $n22s[$i];
  my $x220 = ($n22 > 0)? $n220s[$i]/$n22 : '-';
  my $x221 = ($n22 > 0)? $n221s[$i]/$n22 : '-';

  #   print "$mid  $n00 $x001 $x002  $n01 $x012  $n02 $x020 $x022  $n12 $x120  $n22 $x220 $x221 \n";

  my $n0022 = $n00 + $n22;
  my $n001221 = $n001s[$i] + $n221s[$i];
  my $n002220 = $n002s[$i] + $n220s[$i];
  my $n0112 = $n01 + $n12;
  my $n012120 = $n012s[$i] + $n120s[$i];
  my $n020022 = $n020s[$i] + $n022s[$i];
  
  my $x001221 = ($n0022 > 0)? $n001221/$n0022 : '-';
  my $x002220 = ($n0022 > 0)? $n002220/$n0022 : '-';
  if ($n0022 > 0) {
    push @x001221s, $x001221;
    push @x002220s, $x002220;
  }
  my $x012120 = ($n0112 > 0)? $n012120/$n0112 : '-';
  if ($n0112 > 0) {
    push @x012120s, $x012120;
  }
  my $x020022 = ($n02 > 0)? $n020022/$n02 : '-';
  if ($n02 > 0) {
    push @x020022s, $x020022;
  }
  print "$mid  $n00 $x001221 $x002220  $n0112 $x012120  $n02 $x020022 \n";
}

my @sxs = sort {$a <=> $b} @x001221s;
my $n_prcntl = int ($p * scalar @x001221s);
print "$p  x001221: ", $sxs[$n_prcntl], "\n";

@sxs = sort {$a <=> $b} @x002220s;
$n_prcntl = int ($p * scalar @x002220s);
print "$p  x002220: ", $sxs[$n_prcntl], "\n";

@sxs = sort {$a <=> $b} @x012120s;
$n_prcntl = int ($p * scalar @x012120s);
print "$p  x012120: ", $sxs[$n_prcntl], "\n";

@sxs = sort {$a <=> $b} @x020022s;
$n_prcntl = int ($p * scalar @x020022s);
print "$p  x020022: ", $sxs[$n_prcntl], "\n";


sub good_bad_counts{
  my $acc_gts = shift;		# string
  my $mat_gts = shift;
  my $pat_gts = shift;

  my ($good_count, $bad_count) = (0, 0);
  #    print STDERR "XXXXXXX:   ", length $acc_gts, " ", length $mat_gts, " ", length $pat_gts, "\n";
  for (my $i = 0; $i < length $acc_gts; $i++) {
    my $a_gt = substr($acc_gts, $i, 1);
    my $m_gt = substr($mat_gts, $i, 1);
    my $p_gt = substr($pat_gts, $i, 1);
    if ($m_gt == 1  and  $p_gt == 1) {
      # $a_gt can be 0, 1, or 2. So doesn't add any info. Skip for now.
    } elsif ($m_gt == $p_gt) {	# i.e. either 0,0 or 2,2
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

sub distance{
  my $gstr1 = shift;
  my $gstr2 = shift;
  my $d = 0;
  for (my $i = 0; $i < length $gstr1; $i++) {
    my $g1 = substr($gstr1, $i, 1);
    my $g2 = substr($gstr2, $i, 1);
    $d += abs($g1 - $g2);
  }
  return $d/length $gstr1;
}

sub nxyzs {
  my $acc = shift;
  my $mat = shift;
  my $pat = shift;

  my $L = length $acc;
  die "Length inconsistency \n" if($L != length $mat  or $L != length $pat);
  
  my $n00s = shift;		# array ref.
  my $n001s = shift;
  my $n002s = shift;

  my $n01s = shift;
  my $n012s = shift;

  my $n02s = shift;
  my $n020s = shift;
  my $n022s = shift;

  my $n12s = shift;
  my $n120s = shift;

  my $n22s = shift;
  my $n220s = shift;
  my $n221s = shift;

  my $nX = 0;

  for (my $i = 0; $i < $L; $i++) {
    my $agt = substr $acc, $i, 1;
    my $mgt = substr $mat, $i, 1;
    my $pgt = substr $pat, $i, 1;
    if ($pgt eq '0') {		##### p = 0 #####
      if ($mgt eq '0') {	# 00
	$n00s->[$i]++;
	if ($agt eq '0') {
	} elsif ($agt eq '1') {
	  $n001s->[$i]++;
	} elsif ($agt eq '2') {
	  $n002s->[$i]++;
	} else {
	  $nX++;
	}
      } elsif ($mgt eq '1') {	# 01
	$n01s->[$i]++;
	if ($agt eq '0') {
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	  $n012s->[$i]++;
	} else {
	  $nX++;
	}
      } elsif ($mgt eq '2') {	# 02
	$n02s->[$i]++;
	if ($agt eq '0') {
	  $n020s->[$i]++;
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	  $n022s->[$i]++;
	} else {
	  $nX++;
	}
      } else {
	$nX++;
      }
    } elsif ($pgt eq '1') { 	##### p = 1 #####
      if ($mgt eq '0') {	# 01
	$n01s->[$i]++;
	if ($agt eq '0') {
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	  $n012s->[$i]++;
	} else {
	  $nX++;
	}
      } elsif ($mgt eq '1') {	# 11
	# if ($agt eq '0') {
	#   $n110s->[$i]++;
	# } elsif ($agt eq '1') {
	#   $n111s->[$i]++;
	# } elsif ($agt eq '2') {
	#   $n112s->[$i]++;
	# } else {
	#   $nX++;
	# } 
      } elsif ($mgt eq '2') {	# 12
	$n12s->[$i]++;
	if ($agt eq '0') {
	  $n120s->[$i]++;
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	} else {
	  $nX++;
	}
      } else {
	$nX++;
      }
    } elsif ($pgt eq '2') {	##### p = 2 #####
      if ($mgt eq '0') {	# 02
	$n02s->[$i]++;
	if ($agt eq '0') {
	  $n020s->[$i]++;
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	  $n022s->[$i]++;
	} else {
	  $nX++;
	}
      } elsif ($mgt eq '1') {	# 12
	$n12s->[$i]++;
	if ($agt eq '0') {
	  $n120s->[$i]++;
	} elsif ($agt eq '1') {
	} elsif ($agt eq '2') {
	} else {
	  $nX++;
	}
      } elsif ($mgt eq '2') {	# 22
	$n22s->[$i]++;
	if ($agt eq '0') {
	  $n220s->[$i]++;
	} elsif ($agt eq '1') {
	  $n221s->[$i]++;
	} elsif ($agt eq '2') {
	} else {
	  $nX++;
	}
      } else {
	$nX++;
      }
    } else {
      $nX++;
    }
  }
  return $nX; # ($n00s, $n001s, $n002s,  $n01s, $n012s, $n02s, $n020s, $n022s, $n12s, $n120s, $n22s, $n220s, $n221s, $nX);
}
