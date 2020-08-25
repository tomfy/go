package main

import (
	"bufio"
	//	"container/heap"
	//	"etc"
	"flag"
	"fmt"
	"log"
	"math"
	//	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	//	"strconv"
	"strings"
	"sync"
	//	"time"

	//	"mytypes"
	//	"priorityqueue"
	//	"seqchunkset"
	"sequenceset"
)

type Genotype_diff_counts struct { // store counts of markers with each of 6 possible combinations of genotypes.
	n00     int
	n01     int
	n02     int
	n11     int
	n12     int
	n22     int
	nok     int // sum of n00, n01, n02, n11, n12, n22
	n_other int // count of markers for which one of the gts is something besides 0, 1, 2.
}

type Pedigree_info struct {
	Acc_id       string
	Mat_id       string
	Pat_id       string
	Parental_ids string // ordered, for use as key
	Am_diff      Genotype_diff_counts
	Ap_diff      Genotype_diff_counts
	Mp_diff      Genotype_diff_counts
	Ped_ns       Pedigree_ns
	Ped_xs       Pedigree_xs
}

type Pedigree_ns struct {
	Parental_ids string
	/*	Mat_id       string
		Am_dist      float64
		Pat_id       string
		Ap_dist      float64 */
	n000 int
	n001 int
	n002 int
	n010 int
	n011 int
	n012 int
	n020 int
	n021 int
	n022 int
	n110 int
	n111 int
	n112 int
	n120 int
	n121 int
	n122 int
	n220 int
	n221 int
	n222 int
	nX   int
}

type Pedigree_xs struct {
	Parental_ids string
	/*	acc_id       string
		mat_id       string
		pat_id       string */
	x001 float64
	x002 float64
	x012 float64
	x020 float64
	x022 float64
	x110 float64
	x111 float64
	x112 float64
	x120 float64
	x220 float64
	x221 float64
	/*	am_hgmr      float64
		ap_hgmr      float64 */
}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")
var nfpr = 4 // number of fields per relative (e.g. id agmr hgmr dist  -> nfpr = 4)
// read in
// genotype (matrix) file,
// closest relatives file,
// pedigree file.
// check that parents given in pedigree file are ok, i.e. pedigree_scores for the pedigree is small,
// if not, try various combinations of closest relatives as parents, (n+1 choose 2)
func main() {

	/* command line options: */

	/* input files: */
	var gt_matrix_file, rel_file, ped_file string //
	flag.StringVar(&gt_matrix_file, "gts", "", "filename of matrix-format genotypes file (rows: accessions, columns: markers)")
	flag.StringVar(&rel_file, "rels", "", "filename of close relatives file.")
	flag.StringVar(&ped_file, "peds", "", "filename of pedigree file.")
	

	/* search control parameters */
	var n_close_rels_to_check int
	n_alt_parental_pairs_out := 40
	
	max_md_prop := 0.0
	flag.IntVar(&n_close_rels_to_check, "nrels", 10, "Number of closest distance rels. to check as parents (default: 10)")
	flag.IntVar(&n_alt_parental_pairs_out, "nalts", 40, "Number of alternative parental pairs to output (default: 40)")
	flag.Parse() // parse the command line

	fmt.Fprintf(os.Stderr, "# genotype matrix file: %s pedigree file: %s relatives file: %s\n", gt_matrix_file, ped_file, rel_file)
	// Seed the rng

	// profiling stuff
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	var wg1 sync.WaitGroup
	id_seqset := make(map[string]*sequenceset.Sequence_set) // keeps track of which seq set each seq id belongs to

	//*************** load genotype sequences from file:

	seq_set := &sequenceset.Sequence_set{}
	wg1.Add(1)
	sequenceset.Construct_from_matrix_file(gt_matrix_file, max_md_prop, &id_seqset, seq_set, &wg1)

	n_markers := len(seq_set.Sequences[0])
	fmt.Fprintln(os.Stderr, "# n markers: ", n_markers)
	pedigrees_and_scores := // make(map[string]Pedigree_xs) // pedigrees for pedigree file, and
		make(map[string]*Pedigree_info)

	//****************** read pedigree information from file:
	fh, err := os.Open(ped_file)
	//	fmt.Fprintf(os.Stderr, "# pedigree file name: %s \n", ped_file)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Couldn't open ", ped_file, " for reading.")
		os.Exit(1)
	}
	scanner := bufio.NewScanner(fh)
	scanner.Buffer(make([]byte, 10000), 1000000) // th
	line_number := 0
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line) // split on one or more whitespace chars.
		if line_number == 0 {          // this line should have Accession and other col headings, skip.
			if fields[0] != "Accession" {
				os.Exit(1)
			}
		} else {
			n_fields := len(fields)
			acc_id := fields[n_fields-3]
			mat_id := fields[n_fields-2]
			pat_id := fields[n_fields-1]
			parental_id_string := Ordered_concat_string_pair(mat_id, pat_id)
			if (acc_id != "NA") && (mat_id != "NA") && (pat_id != "NA") { // all three ids must be present, else skip
				//	fmt.Fprintf(os.Stderr, "%s  %s  %s   ", acc_id, mat_id, pat_id)
				acc_gts := seq_set.Sequences[seq_set.SeqId_index[acc_id]]
				mat_gts := seq_set.Sequences[seq_set.SeqId_index[mat_id]]
				pat_gts := seq_set.Sequences[seq_set.SeqId_index[pat_id]]

				ped_ns := calculate_pedigree_ns(parental_id_string, acc_gts, mat_gts, pat_gts)
				ped_xs := calc_ped_xs(parental_id_string, ped_ns)
				pedigrees_and_scores[acc_id] = &Pedigree_info{
					acc_id, mat_id, pat_id, parental_id_string,
					calc_gt_diff_counts(acc_gts, mat_gts),
					calc_gt_diff_counts(acc_gts, pat_gts),
					calc_gt_diff_counts(mat_gts, pat_gts),
					ped_ns,
					ped_xs}
			} else {
				pedigrees_and_scores[acc_id] = &Pedigree_info{Acc_id: acc_id, Mat_id: mat_id, Pat_id: pat_id, Parental_ids: parental_id_string}
			}
		}
		line_number++
	}

	//****************** load closest relatives from file (if specified):
	{
		fh, err := os.Open(rel_file)
		fmt.Fprintf(os.Stderr, "# relatives file name: %s \n", rel_file)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Couldn't open ", rel_file, " for reading.")
			os.Exit(1)
		}
		scanner := bufio.NewScanner(fh)
		scanner.Buffer(make([]byte, 10000), 1000000)

		for scanner.Scan() {
			line := scanner.Text()
			fields := strings.Fields(line) // split on one or more whitespace chars.
			if fields[0] == "#" {
				// comment - skip
			} else {
				n_fields := len(fields)
				acc_id := fields[0]
				if pas, ok := pedigrees_and_scores[acc_id]; ok {

					fmt.Printf("%s  ", pedigree_info_string_x(pas))
					acc_gts := seq_set.Sequences[seq_set.SeqId_index[acc_id]]

					if nfpr*n_close_rels_to_check > n_fields-1 {
						n_close_rels_to_check = (n_fields - 1) / nfpr
					}

					n_parental_pairs_to_try := (n_close_rels_to_check * (n_close_rels_to_check + 1)) / 2
					candidate_pedigrees_info := make([]*Pedigree_info, n_parental_pairs_to_try)
					i := 0
					//	fmt.Fprintf(os.Stderr, "XX: %s \n", line)
				//	fmt.Fprintf(os.Stderr, "AAA:  %d %d \n", n_parental_pairs_to_try, n_close_rels_to_check)
					for k := 0; k < n_close_rels_to_check; k++ {
						mat_id := fields[nfpr*k+1]
						mat_gts := seq_set.Sequences[seq_set.SeqId_index[mat_id]]

						for l := k; l < n_close_rels_to_check; l++ {
							pat_id := fields[nfpr*l+1]
							pat_gts := seq_set.Sequences[seq_set.SeqId_index[pat_id]]
							parental_id_string := Ordered_concat_string_pair(mat_id, pat_id)
							ped_ns := calculate_pedigree_ns(parental_id_string, acc_gts, mat_gts, pat_gts)
							ped_xs := calc_ped_xs(parental_id_string, ped_ns)
							candidate_pedigrees_info[i] = &Pedigree_info{
								acc_id, mat_id, pat_id, parental_id_string,
								calc_gt_diff_counts(acc_gts, mat_gts),
								calc_gt_diff_counts(acc_gts, pat_gts),
								calc_gt_diff_counts(mat_gts, pat_gts),
								ped_ns,
								ped_xs}
							//	fmt.Fprintf(os.Stderr, "%d %d %d %s  ", i, k, l, pedigree_info_string_x(candidate_pedigrees_info[i]))
							i++
						}
					}
					//	fmt.Fprintf(os.Stderr, "\n")

			//		fmt.Fprintln(os.Stderr, "i, len: ", i, len(candidate_pedigrees_info))
					sort.Slice(candidate_pedigrees_info,
						func(i, j int) bool {
							/*	cdii := candidate_pedigrees_info[i].Ped_ns
								cdij := candidate_pedigrees_info[j].Ped_ns
								if true {
									ni := cdii.n020 + cdii.n022
									di := ni + cdii.n021
									nj := cdij.n020 + cdij.n022
									dj := nj + cdij.n021
									return float64(ni)/float64(di) < float64(nj)/float64(dj) // sort from least to greatest
								} else {
									ni := cdii.n002 + cdii.n220
									di := ni + cdii.n000 + cdii.n222
									nj := cdij.n002 + cdij.n220
									dj := nj + cdij.n000 + cdij.n222
									return float64(ni)/float64(di) < float64(nj)/float64(dj) // sort from least to greatest
								} */

							/*	return candidate_pedigrees_info[i].Ped_xs.x002+candidate_pedigrees_info[i].Ped_xs.x220 <
								candidate_pedigrees_info[j].Ped_xs.x002+candidate_pedigrees_info[j].Ped_xs.x220 /* */

							return ped_score_x(candidate_pedigrees_info[i]) < ped_score_x(candidate_pedigrees_info[j])
						})

					naltout := 0
					jped := 200
					_ = jped
					for j, rel_pas := range candidate_pedigrees_info {
						if rel_pas.Parental_ids != pas.Parental_ids {
							fmt.Printf("%s  ", pedigree_info_string_x(rel_pas))
							naltout++
							if naltout >= n_alt_parental_pairs_out {
								break
							}
						} else {
							jped = j
						}
					}
					fmt.Printf("\n")
				//	fmt.Fprintf(os.Stderr, "%s  %s  %d \n", acc_id, pas.Parental_ids, jped)
				} // end of id have pedigree for acc_id
			}
		}
	}
}

func calc_gt_diff_counts(gts1 string, gts2 string) Genotype_diff_counts {
	var n00, n01, n02, n11, n12, n22, nother = 0, 0, 0, 0, 0, 0, 0
	for i := 0; i < len(gts1); i++ {
		c1 := gts1[i : i+1]
		c2 := gts2[i : i+1]
		if c1 == "0" {
			if c2 == "0" {
				n00++
			} else if c2 == "1" {
				n01++
			} else if c2 == "2" {
				n02++
			} else {
				nother++
			}
		} else if c1 == "1" {
			if c2 == "0" {
				n01++
			} else if c2 == "1" {
				n11++
			} else if c2 == "2" {
				n12++
			} else {
				nother++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				n02++
			} else if c2 == "1" {
				n12++
			} else if c2 == "2" {
				n22++
			} else {
				nother++
			}
		} else {
			nother++
		}
	}
	nok := n00 + n01 + n02 + n11 + n12 + n22
	gt_diff_count := Genotype_diff_counts{n00, n01, n02, n11, n12, n22, nok, nother}
	return gt_diff_count
}

// ********************************************************************************
func ahd(gdcs Genotype_diff_counts) (float64, float64, float64) {
	agmr := float64(gdcs.n01+gdcs.n02+gdcs.n12) / float64(gdcs.nok)
	hgmr := float64(gdcs.n02) / float64(gdcs.n00+gdcs.n02+gdcs.n22)
	dist := float64(gdcs.n01+2*gdcs.n02+gdcs.n12) / float64(gdcs.nok)
	return agmr, hgmr, dist
}

func ahd_string(gdcs Genotype_diff_counts) string {
	agmr, hgmr, dist := ahd(gdcs)
	return fmt.Sprintf("%7.5f %7.5f %7.5f", agmr, hgmr, dist)
}

func calc_ped_xs(parents string, pns Pedigree_ns) Pedigree_xs {
	var x001, x002,
		x012,
		x020, x022,
		x110, x111, x112,
		x120,
		x220, x221 float64

	n00 := pns.n000 + pns.n001 + pns.n002
	if n00 > 0 {
		x001 = float64(pns.n001) / float64(n00)
		x002 = float64(pns.n002) / float64(n00)
	} else {
		x001 = 2
		x002 = 2
	}

	n01 := pns.n010 + pns.n011 + pns.n012
	if n01 > 0 {
		x012 = float64(pns.n012) / float64(n01)
	} else {
		x012 = 2
	}

	n02 := pns.n020 + pns.n021 + pns.n022
	if n02 > 0 {
		x020 = float64(pns.n020) / float64(n02)
		x022 = float64(pns.n022) / float64(n02)
	} else {
		x020 = 2
		x022 = 2
	}

	n11 := pns.n110 + pns.n111 + pns.n112
	if n11 > 0 {
		x110 = float64(pns.n110) / float64(n11)
		x111 = float64(pns.n111) / float64(n11)
		x112 = float64(pns.n112) / float64(n11)
	}

	n12 := pns.n120 + pns.n121 + pns.n122
	if n12 > 0 {
		x120 = float64(pns.n120) / float64(n12)
	} else {
		x120 = 2
	}

	n22 := pns.n220 + pns.n221 + pns.n222
	if n22 > 0 {
		x220 = float64(pns.n220) / float64(n22)
		x221 = float64(pns.n221) / float64(n22)
	} else {
		x220 = 2
		x221 = 2
	}
	return Pedigree_xs{
		parents,
		x001, x002,
		x012,
		x020, x022,
		x110, x111, x112,
		x120,
		x220, x221}
}

// returns g_count, b_count pedigree_score = b_count/(g_count + b_count)
func calculate_pedigree_ns(parental_id string, a_gt string, m_gt string, p_gt string) Pedigree_ns { // , int, int, int) {
	nmarkers := len(a_gt)
	if len(m_gt) != nmarkers || len(p_gt) != nmarkers {
		os.Exit(1)
	}
	var g_count, b_count = 0, 0
	var n000, n001, n002,
		n010, n011, n012, // e.g. n010 = number of markers with 0 and 1 in 'parents', 0 in accession.
		n020, n021, n022,
		n110, n111, n112,
		n120, n121, n122,
		n220, n221, n222, nX = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
	for i := 0; i < len(a_gt); i++ {
		a := a_gt[i : i+1]
		m := m_gt[i : i+1]
		p := p_gt[i : i+1]

		if m == "0" {
			if p == "0" {
				if a == "0" {
					g_count++
					n000++
				} else if a == "1" {
					n001++
					b_count++
				} else if a == "2" {
					n002++
					b_count++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "1" {
				if a == "2" {
					b_count++
					n012++
				} else if a == "1" {
					g_count++
					n011++
				} else if a == "0" {
					g_count++
					n010++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "2" {
				if a == "1" {
					g_count++
					n021++
				} else if a == "0" {
					b_count++
					n020++
				} else if a == "2" {
					b_count++
					n022++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else {
				nX++ // p is something besides 0,1,2
			}
		} else if m == "1" {
			if p == "0" {
				if a == "2" {
					b_count++
					n012++
				} else if a == "1" {
					g_count++
					n011++
				} else if a == "0" {
					g_count++
					n010++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "1" {
				//	g_count++ // or don't increment
				if a == "0" {
					n110++
				} else if a == "1" {
					n111++
				} else if a == "2" {
					n112++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "2" {
				if a == "0" {
					b_count++
					n120++
				} else if a == "1" {
					g_count++
					n121++
				} else if a == "2" {
					g_count++
					n122++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else {
				nX++ // p is something besides 0,1,2
			}
		} else if m == "2" {
			if p == "0" {
				if a == "1" {
					g_count++
					n021++
				} else if a == "0" {
					b_count++
					n020++
				} else if a == "2" {
					b_count++
					n022++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "1" {
				if a == "0" {
					b_count++
					n120++
				} else if a == "1" {
					g_count++
					n121++
				} else if a == "2" {
					g_count++
					n122++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else if p == "2" {
				if a == "2" {
					g_count++
					n222++
				} else if a == "0" {
					b_count++
					n220++
				} else if a == "1" {
					b_count++
					n221++
				} else {
					nX++ // a is something besides 0, 1, or 2
				}
			} else {
				nX++ // p is something besides 0,1,2
			}
		} else {
			nX++ // m is something besides 0,1,2
		}
	}

	return Pedigree_ns{parental_id,
		n000, n001, n002,
		n010, n011, n012,
		n020, n021, n022,
		n110, n111, n112,
		n120, n121, n122,
		n220, n221, n222, nX}
}

func pedigree_info_string(pedinf *Pedigree_info) string {
	pxs := pedinf.Ped_xs
	pns := pedinf.Ped_ns
	return fmt.Sprintf("%s  %s %s  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f   %s  %s  %s  %4d %4d %4d %4d %4d %4d %7.5f",
		pedinf.Acc_id,      // col 1 (24)
		pedinf.Mat_id,      // col 2
		pedinf.Pat_id,      // col 3
		pxs.x001, pxs.x002, // 4,5
		pxs.x012,           // 6
		pxs.x020, pxs.x022, // 7,8
		pxs.x110, pxs.x111, pxs.x112, // 9,10,11
		pxs.x120,           // 12
		pxs.x221, pxs.x220, // 13,14
		ahd_string(pedinf.Am_diff), // 15-17  (38-40, ...)
		ahd_string(pedinf.Ap_diff), // cols 18-20
		ahd_string(pedinf.Mp_diff), // 21-23
		pns.n000+pns.n001+pns.n002,
		pns.n010+pns.n011+pns.n012,
		pns.n020+pns.n021+pns.n022,
		pns.n110+pns.n111+pns.n112,
		pns.n120+pns.n121+pns.n122,
		pns.n220+pns.n221+pns.n222,
		ped_score(pedinf),
	)
}

func pedigree_info_string_x(pedinf *Pedigree_info) string {
	//	pxs := pedinf.Ped_xs
	pns := pedinf.Ped_ns

	n00 := pns.n000 + pns.n001 + pns.n002
	n01 := pns.n010 + pns.n011 + pns.n012
	n02 := pns.n020 + pns.n021 + pns.n022
	n11 := pns.n110 + pns.n111 + pns.n112
	n12 := pns.n120 + pns.n121 + pns.n122
	n22 := pns.n220 + pns.n221 + pns.n222

	x012210 := -1.0
	if (n01 + n12) > 1 {
		x012210 = float64(pns.n012+pns.n120) / float64(n01+n12)
	}
	x020022 := -1.0
	if n02 > 0 {
		x020022 = float64(pns.n020+pns.n022) / float64(n02)
	}
	return fmt.Sprintf("%s  %s %s  %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f",
		pedinf.Acc_id, // col 1, 8, 15, ...
		pedinf.Mat_id, // col 2
		pedinf.Pat_id, // col 3
		float64(pns.n001+pns.n221)/float64(n00+n22), // 4
		float64(pns.n002+pns.n220)/float64(n00+n22), // 5
		x012210, // 6
		x020022, // 7
		float64(pns.n111)/float64(n11),
			ped_score_x(pedinf),
	)
}

func ped_score_x(pedinf *Pedigree_info) float64 { // small score suggests parent pair is correct.
	pns := pedinf.Ped_ns
	//	pns := pedinf.Ped_ns
	score := float64(0)
	count := 0
	f := 1.5
	d00 := float64(pns.n000 + pns.n001 + pns.n002 + pns.n220 + pns.n221 + pns.n222)
	d02 := float64(pns.n020 + pns.n021 + pns.n022)
	d12 := float64(pns.n010 + pns.n011 + pns.n012 + pns.n120 + pns.n121 + pns.n122)
	var sink float64 
	if d00 > 0 {
		s002 := float64(pns.n002+pns.n220) / d00
		sink = 0.5 * (1.0 + math.Tanh(f*(s002-0.002)/0.001))
	//	fmt.Fprintf(os.Stderr, "XX  %7.5f ", sink)
		score += sink
		s001 := float64(pns.n001+pns.n221) / d00
		sink := 0.5 * (1.0 + math.Tanh(f*(s001-0.05)/0.01))
		//	fmt.Fprintf(os.Stderr, "%7.5f ", sink)
		score += sink //0.5 * (1.0 + math.Tanh(f*(s001-0.05)/0.015))
		count += 2
	}else{
	//	fmt.Fprintf(os.Stderr, "XX  %7.5f %7.5f ", -1.0, -1.0)
	}
	if d02 > 0 {
		s020 := float64(pns.n020+pns.n022) / d02
	//	fmt.Fprintf(os.Stderr, "xxx %7.5f  \n", s020);
		sink = 0.5 * (1.0 + math.Tanh(f*(s020-0.14)/0.03))
	//		fmt.Fprintf(os.Stderr, "%7.5f ", sink)
		score += sink //0.5 * (1.0 + math.Tanh(f*(s020-0.014)/0.003))
		count++
	}else{
	//	fmt.Fprintf(os.Stderr, "%7.5f ", -1.0)
	}
	if d12 > 0 {
		s012 := float64(pns.n012 + pns.n120) / d12
		sink = 0.5 * (1.0 + math.Tanh(f*(s012-0.12)/0.02))
		//	fmt.Fprintf(os.Stderr, "%7.5f ", sink)
		score += sink
		count++
	}else{
	//	fmt.Fprintf(os.Stderr, "%7.5f ", -1.0)
	}
//	fmt.Fprintf(os.Stderr, "\n")
	return score/float64(count)
}

func ped_score(pedinf *Pedigree_info) float64 {
	pxs := pedinf.Ped_xs
	//	pns := pedinf.Ped_ns
	score := float64(0)
	count := 0
	f := 1.5
	if pxs.x001 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x001-0.05)/0.01))
		count++
	}
	if pxs.x002 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x002-0.002)/0.001))
		count++
	}
	if pxs.x012 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x012-0.02)/0.007))
		count++
	}
	if pxs.x020 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x020-0.07)/0.01))
		count++
	}
	if pxs.x022 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x022-0.06)/0.01))
		count++
	}
	if pxs.x120 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x120-0.025)/0.007))
		count++
	}
	if pxs.x220 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x220-0.02)/0.01))
		count++
	}
	if pxs.x221 < 1.0 {
		score += 0.5 * (1.0 - math.Tanh(f*(pxs.x221-0.08)/0.02))
		count++
	}
	if count > 0 {
		score /= float64(count)
	}
	return 1.0 - score //, count
}

// get the query id and seq, and the candidate match ids and seqs from channel, calculate distances
// send them through a channel
/* func calculate_distances(inch chan []*mytypes.IdSeq, outch chan []*mytypes.IdCmfDistance, wg *sync.WaitGroup) {
	defer wg.Done()
	for {
		q_and_matches, ok := <-inch
		if !ok {
			break
		}
		size := len(q_and_matches)
		out := make([]*mytypes.IdCmfDistance, size)
		query := q_and_matches[0]
		qout := mytypes.IdCmfDistance{query.Id, -1, -1} // set Distance negative here to guarantee sort puts query in position 0
		out[0] = &qout
		for i := 1; i < len(q_and_matches); i++ {
			subj := q_and_matches[i]

			sout := mytypes.IdCmfDistance{subj.Id, -1, dist}
			out[i] = &sout
		}
		sort.Slice(out, func(i, j int) bool { // sort them by ChunkMatchFraction
			return out[i].Distance < out[j].Distance
		})
		out[0].Distance = 0 // now set it to 0
		outch <- out
	}
} /* */

/*  func  calculate_candidate_distances(id_seqset map[string]*sequenceset.Sequence_set, qid_cmfpq map[string]*priorityqueue.PriorityQueue) int {
	dist_calc_count := 0
	idpair_dist := make(map[string]float64)
	for id1, cmfpq := range qid_cmfpq {
		seqset1 := id_seqset[id1]
		seq1 := seqset1.Sequences[seqset1.SeqId_index[id1]]
		top_matches := make([]mytypes.IdCmfDistance, len(*cmfpq))
		fmt.Printf("%s ", id1)
		for i, cmf := range *cmfpq {
			id2 := cmf.Id
			cmf := cmf.Cmf
			seqset2 := id_seqset[id2]
			seq2 := seqset2.Sequences[seqset2.SeqId_index[id2]]

			var idpair string
			if id1 < id2 { // put ids together with the 'smaller' one on left:
				idpair = id1 + "\t" + id2
			} else {
				idpair = id2 + "\t" + id1
			}
			dist, ok := idpair_dist[idpair]
			if !ok {                        // don't already have the distance for this pair - calculate it from the seq1, seq2.
				n00_22, n11, nd1, nd2 := sequenceset.Distance(seq1, seq2)
				dist = float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
				idpair_dist[idpair] = dist
				dist_calc_count++
			}
			icd := mytypes.IdCmfDistance{id2, cmf, dist}
			top_matches[i] = icd
		}
		sort.Slice(top_matches, func(i, j int) bool {
			return top_matches[i].Distance < top_matches[j].Distance
		})
		for _, x := range top_matches {
			//	fmt.Printf("%s %5.3f %7.5f  ", x.Id, x.ChunkMatchFraction, x.Distance)
			fmt.Printf("%s %7.5f  ", x.Id, x.Distance)
		}
		fmt.Println()
	}
	return dist_calc_count
} /* */

func what_is_format(filename string) string { // returns "matrix", "fasta" or "other"
	//	fmt.Println("filename: ", filename)
	fh, err := os.Open(filename)
	if err != nil {
		fmt.Println("Couldn't open ", filename)
		os.Exit(1)
	}
	scanner := bufio.NewScanner(fh)
	//	buf := make([]byte, 10000)
	scanner.Buffer(make([]byte, 10000), 1000000) // the default here was 64*1024 - not big enough! Prob. need to
	//	fmt.Printf("max token size of scanner: %d\n", scanner.maxTokenSize)
	//	fmt.Printf("scanner.Scan() returned: %t \n", scanner.Scan())
	//	fmt.Println("[" + scanner.Text() +"]")
	for scanner.Scan() {
		line := scanner.Text()
		//	fmt.Println("Line: [" + line[0:50])
		fields := strings.Fields(line)
		//	fmt.Println("first line, first field: ", fields[0])
		if fields[0] == "MARKER" {
			return "matrix"
		} else if fields[0][0:1] == ">" {
			return "fasta"
		} else {
			fmt.Println("Unknown file format. ", fields[0])
			return "other"
		}
		break
	}
	return "other"
}

func Ordered_concat_string_pair(a string, b string) string { // returns concatenation of a and b ordered such that the 'lesser' one is first.
	if strings.Compare(a, b) != 1 {
		return a + " " + b
	} else {
		return b + " " + a
	}
}

/* func Initialize_priorityqueues(seq_set *sequenceset.Sequence_set, qid_cmfpq *map[string]*priorityqueue.PriorityQueue) {
	// make an empty priority queue for each sequence to hold the best candidate matches to that sequence
	for qid, _ := range seq_set.SeqId_index {
		pq := make(priorityqueue.PriorityQueue, 0)
		heap.Init(&pq)
		(*qid_cmfpq)[qid] = &pq
	}
} /* */

func WaitForEnter(message string, actually_wait bool) {
	//	fmt.Fprintln(os.Stderr, "warning msg: ", message)
	//	fmt.Fprintln(os.Stderr, "actually_wait:  ", actually_wait)
	if actually_wait {
		reader := bufio.NewReader(os.Stdin)
		fmt.Fprintln(os.Stderr, message)
		x, _ := reader.ReadString('\n')
		_ = x
	}
}

func calculate_agmr_hgmr(seq1 string, seq2 string) (float64, float64) {
	n00_22, n11, nd1, nd2 := sequenceset.Distance(seq1, seq2)
	agmr := float64(nd1+nd2) / float64(n00_22+n11+nd1+nd2)
	hgmr := float64(nd2) / float64(n00_22+nd2)
	return agmr, hgmr
}

func MemUsageString() string {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	result := ""
	// For info on each, see: https://golang.org/pkg/runtime/#MemStats
	result += fmt.Sprintf("# Allocated heap objects: %v MiB;  ", (m.Alloc)/1024/1024)
	//   fmt.Printf("\tTotalAlloc = %v MiB", (m.TotalAlloc)/1024/1024)
	result += fmt.Sprintf("Sys = %v MiB;  ", (m.Sys)/1024/1024)
	result += fmt.Sprintf("Number of garbage collections: %v;  ", m.NumGC)
	result += fmt.Sprintf("Mallocs: %v  Frees: %v  Mallocs-Frees: %v ", m.Mallocs, m.Frees, m.Mallocs-m.Frees)
	return result
}

/* func store_matches(ch chan map[string][]*mytypes.MatchInfo, qid_allokmatches map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	for {
		qid_okmatches, ok := <-ch // ok will be false iff channel is empty and closed.
		if ok {
			for qid, okmatches := range qid_okmatches {
				//	fmt.Fprintln(os.Stderr, "type of: ", reflect.TypeOf(okmatches))
				x, ok := qid_allokmatches[qid]
				if ok {
					qid_allokmatches[qid] = append(x, okmatches...)
				} else {
					qid_allokmatches[qid] = okmatches
				}
			}
		} else {
			return
		}
	}

} /* */

// get the top n_keep matches to each q seq
/* func search(qscs *seqchunkset.Sequence_chunk_set, scs *seqchunkset.Sequence_chunk_set, n_keep int, ch chan map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	qid_okmatches, qid_badmatches, total_chunk_match_count, total_mdmd_match_count :=
		scs.Search_qs(qscs, n_keep) //, &qid_cmfpq)
	_ = total_chunk_match_count
	_ = total_mdmd_match_count
	ch <- qid_okmatches
	//	fmt.Fprintln(os.Stderr, "len(qid_okmatches): ", len(qid_okmatches))
	if len(qid_badmatches) > 0 {
		fmt.Println("#  there are this many bad matches: ", len(qid_badmatches))
	}
} /* */

/* func send_top_candidates(qid_matchinfos map[string][]*mytypes.MatchInfo, id_seqset map[string]*sequenceset.Sequence_set, n_keep int, wg *sync.WaitGroup,
	ch chan []*mytypes.IdSeq) {
	defer wg.Done()
	for q_id, okmatches := range qid_matchinfos { // for each query there are n_cpus*n_keep candidates
		sort.Slice(okmatches, func(i, j int) bool { // sort them by ChunkMatchFraction
			return okmatches[i].ChunkMatchFraction > okmatches[j].ChunkMatchFraction
		})

		qss := id_seqset[q_id]
		//	_ = qss.Check_seq_index_id_maps()
		q_seq := qss.Sequences[qss.SeqId_index[q_id]]
		//	top_matches := make([]mytypes.IdCmfDistance, n_keep)

		q_idseq := mytypes.IdSeq{q_id, q_seq}

		n_to_do := etc.MinInt(n_keep, len(okmatches))


		ids_and_seqs := make([]*mytypes.IdSeq, 0, n_to_do+1) // ids_and_seqs[0] will be the query, others the candidate matches
		ids_and_seqs = append(ids_and_seqs, &q_idseq)

		for j := 0; j < n_to_do; j++ { // calculate distance for the n_keep top candidates.
			matchinfo := okmatches[j]
			s_id := matchinfo.Id
			sss := id_seqset[s_id] // get the relevant Sequence_set
			//	_ = sss.Check_seq_index_id_maps()
			//	sidx1 := sss.SeqId_index[s_id]
			s_index := matchinfo.Index
			//	sid2 := sss.SeqIndex_id[sidx2]
			//	fmt.Println("  s id, index1,2, sid2: ", s_id, sidx1, sidx2, sid2)
			s_seq := sss.Sequences[s_index] // sss.Sequences[sss.SeqId_index[s_id]]
			s_idseq := mytypes.IdSeq{s_id, s_seq}
			ids_and_seqs = append(ids_and_seqs, &s_idseq)
		}
		ch <- ids_and_seqs
	} // end loop over queries
}

func output(ch chan []*mytypes.IdCmfDistance) {

	for q_and_matches := range ch {
		query := q_and_matches[0]
		fmt.Printf("%s  ", query.Id)
		for i := 1; i < len(q_and_matches); i++ {
			subj := q_and_matches[i]
			fmt.Printf("%s %7.5f  ", subj.Id, subj.Distance)
		}
		fmt.Println("")
	}
} /* */

// **********************************************************************************************************************

/* func store_matches_x(ch chan sequenceset.QsetSsetQSmi, qid_matchinfos map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	for {
		qsqsmi, ok := <-ch // ok will be false iff channel is empty and closed.
		if ok {
			qset := qsqsmi.Qss
			sset := qsqsmi.Sss
			for qi, s_mi := range qsqsmi.Qs_mi { // s_mi slice of pointers to MatchInfo
				//	fmt.Fprintln(os.Stderr, "type of: ", reflect.TypeOf(okmatches))

				qid := qset.SeqIndex_id[qi]

				minfos, ok := qid_matchinfos[qid]
				if ok { // we already have some matches to this query id
					//	new_minfos := make([]mytypes.MatchInfo, len(s_mi))
					for si, minfo := range s_mi {
						sid := sset.SeqIndex_id[si]
						minfo.Id = sid
						//	new_minfos[si] = minfo
						qid_matchinfos[qid] = append(minfos, minfo)
					}
					//	qid_matchinfos[qid] = append(minfos, new_minfos)
				} else { // no matches yet to this query id, because this is the first set of subjs done.
					for si, minfo := range s_mi {
						sid := sset.SeqIndex_id[si]
						minfo.Id = sid
						//	new_minfos[si] = minfo
					}
					qid_matchinfos[qid] = s_mi

				}
			}
		} else { // channel is empty and closed.
			return
		}
	}
	return
}
func search_x(qscs *seqchunkset.Sequence_chunk_set, sscs *seqchunkset.Sequence_chunk_set, n_keep int, ch chan sequenceset.QsetSsetQSmi, wg *sync.WaitGroup) {
	defer wg.Done()
	qs_matchinfos := sscs.Search_qs_x(qscs, n_keep)
	x := sequenceset.QsetSsetQSmi{qscs.Sequence_set, sscs.Sequence_set, qs_matchinfos}
	ch <- x
} /* */

/* func distance_old(seq1 string, seq2 string) float64 {
	zero_count := 0
	one_count := 0
	two_count := 0
	for i := 0; i < len(seq1); i++ {
		c1 := seq1[i : i+1]
		c2 := seq2[i : i+1]
		if c1 == "0" {
			if c2 == "0" {
				zero_count++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				two_count++
			}
		} else if c1 == "1" {
			if c2 == "0" {
				one_count++
			} else if c2 == "1" {
				zero_count++
			} else if c2 == "2" {
				one_count++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				two_count++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				zero_count++
			}
		}
	}
	ok_count := zero_count + one_count + two_count // number of sites where neither seq has missing data
	dist_count := one_count + 2*two_count          // sums differences, i.e. 0-1 -> +=1, 0-2 -> += 2, ...
	var distance float64
	if ok_count > 0 {
		distance = float64(dist_count) / float64(ok_count)
	} else {
		distance = -1.0 // couldn't calculate because no sites without missing data
	}
	return distance
} /* */

/* func distance_z(seq1 string, seq2 string) (int, int, int, int) {
	//	zero_count := 0
	one_count := 0
	two_count := 0
	n00_22 := 0 // homozygous, no change, i.e. 0->0 or 2->2
	n11 := 0    // heterozygous, no change, i.e. 1 -> 1
	// n02 := 0 // same as two_count
	// n01_12 := 0 // same as one_count
	for i := 0; i < len(seq1); i++ {
		c1 := seq1[i : i+1]
		c2 := seq2[i : i+1]
		if c1 == "0" {
			if c2 == "0" {
				//	zero_count++
				n00_22++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				two_count++
			}
		} else if c1 == "1" {
			if c2 == "0" {
				one_count++
			} else if c2 == "1" {
				//	zero_count++
				n11++
			} else if c2 == "2" {
				one_count++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				two_count++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				//	zero_count++
				n00_22++
			}
		}
	}
	return n00_22, n11, one_count, two_count
} /* */

// version using sequences stored as slices, rather than as strings.
/* func distance_x(seq1 []uint, seq2 []uint) float64 {
	ok_count := 0   // counts sites where neither seq has missing data
	dist_count := 0 // sums differences, i.e. 0-1 -> +=1, 0-2 -> += 2, ...
	for i := 0; i < len(seq1); i++ {
		c1 := seq1[i]
		c2 := seq2[i]
		if c1 == 0 {
			if c2 == 0 {
				ok_count++
			} else if c2 == 1 {
				ok_count++
				dist_count++
			} else if c2 == 2 {
				ok_count++
				dist_count += 2
			}
		} else if c1 == 1 {
			if c2 == 0 {
				ok_count++
				dist_count++
			} else if c2 == 1 {
				ok_count++
			} else if c2 == 2 {
				ok_count++
				dist_count++
			}
		} else if c1 == 2 {
			if c2 == 0 {
				ok_count++
				dist_count += 2
			} else if c2 == 1 {
				ok_count++
				dist_count++
			} else if c2 == 2 {
				ok_count++
			}
		}
	}
	var distance float64
	if ok_count > 0 {
		distance = float64(dist_count) / float64(ok_count)
	} else {
		distance = -1.0 // couldn't calculate because no sites without missing data
	}
	return distance
}
*/
