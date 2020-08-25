#!/usr/bin/perl -w
use strict;

# read in genotypes file
# rows: accessions, columns: markers
# and pedigreetest output file
# with pedigree id-triple (acc_id parent1_id parent2_id)
# and id-triples for several alternative possible parentages
# and output colony input format file for each triple

my $peditest_filename = shift;
my $accs_to_do = shift // 1;
my $ptfile_span = 9; # the number of columns for each parentage hypothesis.


my %id_gtstr = ();
my %markerid_012count = ();
my $first_line = <>;
my @cols = split(" ", $first_line);
die "first line doesnt have MARKER at left.\n" if ($cols[0] ne 'MARKER');
shift @cols;
my $n_markers = scalar @cols;
my @marker_ids = @cols;

my $accession_number = 0;
my %alphaid_numid = ();		# 
my %numid_alphaid = ();
while (<>) {
  @cols = split(" ", $_);
  my $accid = shift @cols;
  $alphaid_numid{$accid} = $accession_number;
  $numid_alphaid{$accession_number} = $accid;
  $accession_number++;
  while (my ($i, $gt) = each @cols) {
    my $mid = $marker_ids[$i];
    if (!exists $markerid_012count{$mid}) {
      $markerid_012count{$mid} = [0, 0, 0];
    }
    $markerid_012count{$mid}->[$gt]++;
  }
  my $gtstr = join('', @cols);
  $id_gtstr{$accid} = $gtstr;
}


my $h_string = '';
$h_string .= "'colony_temp_out'  ! Dataset name\n";
$h_string .= "'colony_temp_out'  ! Output file name\n";
$h_string .= "1           ! Number of offspring in the sample\n";
$h_string .= "$n_markers        ! Number of loci\n";
$h_string .= "1234        ! RNG seed\n";
$h_string .= "0           ! 0/1=Not updating / update allele frequency\n";
$h_string .= "2           ! 2/1= Dioecious/Monoecious species\n";
$h_string .= "1           ! 0/1=No inbreeding/inbreeding\n";
$h_string .= "0           ! 0/1=Diploid species/HaploDiploid species\n";
$h_string .= "0  0        ! 0/1=Polygamy/Monogamy for males & females\n";
$h_string .= "0           ! 0/1=Clone inference =No/Yes\n";
$h_string .= "1           ! 0/1=Full sibship size scaling =No/Yes\n"; # ???
$h_string .= "1 1.0 1.0   ! 0,1,2,3=No,weak,medium,strong sibship size prior; mean paternal & meteral sibship size\n";
$h_string .= "1           ! 0/1=Unknown/Known population allele frequency\n"; # perhaps better to have '1'
# 12 13 14 15 16  !Number of alleles per locus
my $n_allele_string = '';
my $allele_freq_str = '';
while (my($i, $mid) = each @marker_ids) {
  $n_allele_string .= "2 ";
  my $af1  = allele_freqs( $markerid_012count{$mid} );
  my $af2 = 1.0 - $af1;
  $allele_freq_str .= "1 2\n" . "$af1 $af2\n";
}
$h_string .= "$n_allele_string  ! Number of alleles per locus\n";
$h_string .= $allele_freq_str;

$h_string .= "\n";

$h_string .= "1         ! Number of runs\n";
$h_string .= "2         ! Length of run (1: short, 2: medium, etc.)\n";
$h_string .= "0         ! 0/1=Monitor method by Iterate#/Time in second\n";
$h_string .= "100000    ! Monitor interval in Iterate# / in seconds\n";
$h_string .= "0         ! non-Windows version\n";
$h_string .= "1         ! Full likelihood\n";
$h_string .= "3         ! 1/2/3=low/medium/high Precision for Full likelihood\n";
$h_string .= "mrkr@     ! generic marker name (names will be mrkr1, mrkr2, etc.\n";
$h_string .= "0@        ! all markers codominant\n";
$h_string .= "0.01@     ! allelic dropout rate (for all markers)\n";
$h_string .= "0.01@     ! other genotyping error rate (for all markers)\n";

open my $fh, "<", "$peditest_filename" or die;
my $i = 0;
while (my $line = <$fh>) {
  my $output_filename = "colony_input_" . $i;
  my ($a_gtstr, $p_gtstr, $z_str) = colony_output_string($line, \%id_gtstr);
  $i++;
  $h_string .= "$a_gtstr" . "\n";
  $h_string .= "0.5 0.5 ! prob. of dad/mum included in candidates.\n";
  $h_string .= "\n";

  $h_string .= "$p_gtstr";
  $h_string .= $z_str;
  last if($i >= $accs_to_do);
}
# exit;


open my $fh1, ">", "colony_temp_in";
print $fh1 ($h_string, "\n");

my $colony_stdout = `colony2s IFN:colony_temp_in`; # outputs files colony_temp_out.
my $parentage_filename = "colony_temp_out.ParentPair";

open my $fhz, "<", "$parentage_filename";
my $line1 = <$fhz>;
my %amp_p = ();
if ( $line1 =~ /^\s*OffspringID,InferredDad,InferredMum,Probability\s*$/ ) {
  while (my $line = <$fhz>) {
    $line =~ s/\s+//g;
    my ($aid, $mid, $pid, $prob) = split(',', $line);
 #   print STDERR "$line\n";
 #   print STDERR "$aid $mid $pid $prob\n";
    $aid =~ s/[O]_//;
    $mid =~ s/[MP]_//;
    $pid =~ s/[MP]_//;
#    print STDERR "$aid $mid $pid $prob\n";
    $aid = $numid_alphaid{$aid};
    $mid = $numid_alphaid{$mid};
    $pid = $numid_alphaid{$pid};
 #   print STDERR "$aid $mid $pid $prob\n";
    my $mpids = (($mid cmp $pid) <= 0)? "$mid $pid" : "$pid $mid";
    my $idtriple = "$aid $mpids";
#    print STDERR "idtriple: $idtriple \n";
    $amp_p{$idtriple} += $prob;
  }
}
# print STDERR "n triples: ", scalar keys %amp_p, "\n";

my @sorted_idtriples = sort {$amp_p{$b} <=> $amp_p{$a}} keys %amp_p;
for my $idtriple (@sorted_idtriples) {
  print "$idtriple  ", $amp_p{$idtriple}, "\n";
}



sub colony_output_string{
  my $line = shift;
  my $id_gtstr = shift;
  my @cols = split(" ", $line);
  my $aid = $cols[0];
  my $acc_numid = $alphaid_numid{$aid} // die;
#  print STDERR "aid: $aid  $acc_numid\n";
  my ($a_size, $a_gtstr) = gts_to_allele_pairs( $id_gtstr->{$aid} );
  #  my $allelepairs_a = gts_to_allele_pairs($a_gtstr);
  my $a_str = "O_" . "$acc_numid  " . $a_gtstr . "\n";
  my %parental_ids = ();
  my $parentpair_candidate_count = 0; # including pedigree
  while (@cols) {
    my ($aid1, $mid, $pid) = @cols[0,1,2];
    #    print STDERR "aid1: $aid1 mid: $mid pid: $pid \n";
    die "accession id inconsistency:  $aid $aid1 \n" if($aid1 ne $aid);
    $parental_ids{$mid} = $id_gtstr->{$mid};
    $parental_ids{$pid} = $id_gtstr->{$pid};
    # @cols = @cols[$ptfile_span:-1];
    splice @cols, 0, $ptfile_span;
    $parentpair_candidate_count++;
  }
  my $n_parent_candidates = scalar keys %parental_ids;
  print "# n candidate parent pairs: $parentpair_candidate_count   n distinct indiv. parent candidates: $n_parent_candidates. \n";
  my $parental_str = "$n_parent_candidates $n_parent_candidates  ! numbers of candidates, males & females\n\n"; #
  my ($m_str, $p_str) = ('', '');
  my $id_number = 0;
  for my $id (keys %parental_ids) {
    my $gtstr = $id_gtstr->{$id};
    $id_number = $alphaid_numid{$id} // die "Id $id unknown??\n";
    my ($p_size, $p_gtstr) =  gts_to_allele_pairs($gtstr);
    $m_str .= "M_" . "$id_number  $p_gtstr\n";
    $p_str .= "P_" . "$id_number  $p_gtstr\n";
    #  $id_number++;
  }
  $parental_str .= $m_str . "\n" . $p_str;

  my $z_str =  "0  0          !#known mother-offspring dyads, maternity exclusion threshold\n";
  $z_str .= "0  0          !#known father-offspring dyads, paternity exclusion threshold\n";
  $z_str .= "0             !#known paternal sibship with unknown fathers\n";
  $z_str .= "0             !#known maternal sibship with unknown mothers\n";
  $z_str .= "0             !#known paternity exclusions\n";
  $z_str .= "0             !#known maternity exclusions\n";
  $z_str .= "0             !#known paternal sibship exclusions\n";
  $z_str .= "0             !#known maternal sibship exclusions\n";

  return ($a_str, $parental_str, $z_str);
}

sub gts_to_allele_pairs{
  my $gtstr = shift;
  my @gts = split('', $gtstr);
  my $outstr = '';
  for my $agt (@gts) {
    if ($agt eq '0') {
      $outstr .= '1 1  ';

    } elsif ($agt eq '1') {
      $outstr .= '1 2  ';
    } elsif ($agt eq '2') {
      $outstr .= '2 2  ';
    } else {
      die "Unknown genotype: $agt \n";
    }
  }
  return (scalar @gts, $outstr);
}

sub allele_freqs{
  my $gt_count = shift;		# e.g. [89,64,20]
  my ($n0, $n1, $n2) = @$gt_count;
  my $allele1_count = 2*$n0 + $n1;
  my $allele2_count = $n1 + 2*$n2;
  my $allele1_freq = $allele1_count/($allele1_count + $allele2_count);
  return $allele1_freq;
}
