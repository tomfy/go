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
my @ids = keys %id_genotypes;
my $n_ids = scalar @ids;

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
      my $rid1 = $ids[rand($n_ids)];
      my $r1_gts = $id_genotypes{$rid1};
      my $r2_gts = $id_genotypes{$ids[rand($n_ids)]};
      if (defined $a_gts  and  defined $m_gts  and  defined $p_gts) {
         # score rel. to parents specified in the pedigree.
         my ($goods, $bads) = good_bad_counts($a_gts, $m_gts, $p_gts);
         my $score = $bads/($goods + $bads);
         print STDERR "Score: $score \n";

         my $am_dist = distance($a_gts, $m_gts);
         my $ap_dist = distance($a_gts, $p_gts);
         my $ar1_dist = distance($a_gts, $r1_gts);

         my $ok_parent_id = ($ap_dist < $am_dist)? $pat_id : $mat_id;
          
         my $best_score = 100000;
         my $best_id = 'X';
         if ($score > 0.3) {   # i.e. a 'bad' score
             print STDERR "ok_parent_id: $ok_parent_id \n";
            for my $this_id (@ids) {
               if ($this_id ne $acc_id) {
                  my $this_gts = $id_genotypes{$this_id};
                  my ($gs, $bs) = good_bad_counts($a_gts, $id_genotypes{$ok_parent_id}, $this_gts);
                  my $this_score = $bs/($gs + $bs);
                  if ($this_score < $best_score) {
                     $best_score = $this_score;
                     $best_id = $this_id;
                  }
               }
            }
         }

         my ($hr_goods, $hr_bads) = good_bad_counts($a_gts, $m_gts, $r1_gts);
         my $hr_score = $hr_bads/($hr_goods + $hr_bads);

         my ($r_goods, $r_bads) = good_bad_counts($a_gts, $r1_gts, $r2_gts);
         my $r_score = $r_bads/($r_goods + $r_bads);

         print STDERR "$acc_id  $mat_id  $pat_id    $goods  $bads   $score   $hr_score $r_score  $am_dist $ap_dist $ar1_dist $best_score  $ok_parent_id $best_id \n";
         print "$acc_id  $mat_id  $pat_id    $goods  $bads   $score   $hr_score $r_score  $am_dist $ap_dist $ar1_dist $best_score  $ok_parent_id $best_id \n";
   
      } else {
         warn "One of the 3 ids (accession of parents) not present in genotypes matrix file.\n"
      }
   }


}

sub good_bad_counts{
   my $acc_gts = shift;         # string
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
      } elsif ($m_gt == $p_gt) { # i.e. either 0,0 or 2,2
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
      } else {                  # $p_gt == 2
         return ($a_gt == 1)? 1 : 0;
      }
   } elsif ($m_gt == 1) {
      if ($p_gt == 0) {
         return ($a_gt == 2)? 0 : 1;
      } elsif ($p_gt == 1) {
         return 1;
      } else {                  # $p_gt == 2
         return ($a_gt == 0)? 0 : 1;
      }
      
   } else {			# $m_gt == 2
      if ($p_gt == 0) {
         return ($a_gt == 1)? 1 : 0;
      } elsif ($p_gt == 1) {
         return ($a_gt == 0)? 0 : 1;
      } else {                  # $p_gt == 2
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
