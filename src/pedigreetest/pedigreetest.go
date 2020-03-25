package main

import (
	"bufio"
	"container/heap"
	"etc"
	"flag"
	"fmt"
	"log"
//	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strconv"
	"strings"
	"sync"
//	"time"

	"mytypes"
	"priorityqueue"
	"seqchunkset"
	"sequenceset"
)

type Pedigree_info struct {
	Mat_id             string
	Am_dist            float64
	Pat_id             string
	Ap_dist            float64
	Amp_pedigree_score float64
}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

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
	var n_close_rels_to_check  int
	//	var pauses bool
//	var seed int64
	max_md_prop := 0.0
//	flag.IntVar(&n_cpus, "cpus", 1, "Number of cpus to use.")
//	flag.IntVar(&n_keep, "keep", 20, "Number of best matches to keep.")
//	flag.Int64Var(&seed, "seed", -1, "Rng seed (default: set from clock.)")
	flag.IntVar(&n_close_rels_to_check, "nrels", 5, "Number of closest distance rels. to check as parents (default: 5)")
	/*	flag.IntVar(&n_reps, "reps", 1, "Number of times to repeat the whole search with different random chunk sets.")
		var missing_data_prob, max_missing_data_proportion float64
		flag.Float64Var(&missing_data_prob, "miss", -1, "Proportion missing data to add to genotypes.")
		flag.Float64Var(&max_missing_data_proportion, "max_md", 0.1, "Max proportion of missing data to use marker in chunk set.")
		flag.BoolVar(&pauses, "pause", false, "false: no pauses, true: pause at various points until enter key pressed.") /* */
	//	max_missing_data_proportion := 0
	flag.Parse() // parse the command line

	fmt.Fprintf(os.Stderr, "%s %s %s\n", gt_matrix_file, ped_file, rel_file)
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

	//	var seq_set *sequenceset.Sequence_set
	seq_set := &sequenceset.Sequence_set{}
	//	sequenceset.Construct_sets_from_matrix_file(gt_matrix_file, n_cpus, max_missing_data_proportion, &id_seqset, seq_set)
	wg1.Add(1)
	sequenceset.Construct_from_matrix_file(gt_matrix_file, max_md_prop, &id_seqset, seq_set, &wg1)

	n_markers := len(seq_set.Sequences[0])
	fmt.Fprintln(os.Stderr, "# n markers: ", n_markers)

	// make(map[string]int)
	//	amp_infos := make(map[string]*Pedigree_infoAmp_info,

	pedigrees_and_scores := make(map[string]Pedigree_info)
	//****************** read pedigree information from file:
	fh, err := os.Open(ped_file)
	fmt.Fprintf(os.Stderr, "pedigree file name: %s \n", ped_file)
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
			if acc_id != "NA" && mat_id != "NA" && pat_id != "NA" { // all three ids must be present, else skip
				fmt.Fprintf(os.Stderr, "%s  %s  %s\n", acc_id, mat_id, pat_id)
				acc_gt := seq_set.Sequences[seq_set.SeqId_index[acc_id]]
				mat_gt := seq_set.Sequences[seq_set.SeqId_index[mat_id]]
				pat_gt := seq_set.Sequences[seq_set.SeqId_index[pat_id]]
				g_count, b_count := calculate_pedigree_score(acc_gt, mat_gt, pat_gt)
				fmt.Fprintf(os.Stderr, "%8d %8d \n", g_count, b_count)
				pedigree_score := float64(b_count) / float64(g_count+b_count)
				//	fmt.Println(acc_id, mat_id, pat_id, pedigree_score)
				n00_22, n11, nd1, nd2 := sequenceset.Distance(acc_gt, mat_gt)
				am_dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
				n00_22, n11, nd1, nd2 = sequenceset.Distance(acc_gt, pat_gt)
				ap_dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
				pedigrees_and_scores[acc_id] = Pedigree_info{mat_id, am_dist, pat_id, ap_dist, pedigree_score}
			}
		}
		line_number++
	}

	/* for _, pas := range pedigrees_and_scores {
		fmt.Println("X: ", acc_id, pas[acc_id].Mat_id, pas[acc_id].Pat_id, pas[acc_id].Pedigree_score)
	} /* */

	//****************** load closest relatives from file (if specified):
	{
		fh, err := os.Open(rel_file)
		fmt.Fprintf(os.Stderr, "relatives file name: %s \n", rel_file)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Couldn't open ", rel_file, " for reading.")
			os.Exit(1)
		}
		scanner := bufio.NewScanner(fh)
		scanner.Buffer(make([]byte, 10000), 1000000) // th
	//	line_number := 0
		for scanner.Scan() {
			line := scanner.Text()
			fields := strings.Fields(line) // split on one or more whitespace chars.
			if fields[0] == "#" {
				// comment - skip
			} else {
				n_fields := len(fields)
				acc_id := fields[0]
			//	fmt.Println("Accession: ", acc_id)
				acc_gt := seq_set.Sequences[seq_set.SeqId_index[acc_id]]

			//	n_close_rels_to_check := 7
				if 2*n_close_rels_to_check > n_fields-1 {
					n_close_rels_to_check = (n_fields - 1) / 2
				}
				/*		close_rel_ids := make([]string, n_close_rels_to_check)
						for j := 0; j < n_close_rels_to_check; j++ {
							close_rel_ids[j] = fields[2*j+1]
						} /* */
				n_parental_pairs_to_try := (n_close_rels_to_check * (n_close_rels_to_check + 1)) / 2
				candidate_pedigrees_info := make([]Pedigree_info, n_parental_pairs_to_try)
				i := 0
				for k := 0; k < n_close_rels_to_check; k++ {
					mat_id := fields[2*k+1]
					mat_gt := seq_set.Sequences[seq_set.SeqId_index[mat_id]]

					am_dist, _ := strconv.ParseFloat(fields[2*k+2], 64)
					for l := k; l < n_close_rels_to_check; l++ {
						pat_id := fields[2*l+1]
						pat_gt := seq_set.Sequences[seq_set.SeqId_index[pat_id]]
						ap_dist, _ := strconv.ParseFloat(fields[2*l+2], 64) //fields[2*l+2]
						good_count, bad_count := calculate_pedigree_score(acc_gt, mat_gt, pat_gt)
						p_score := float64(bad_count) / float64(good_count+bad_count)
						p_info := Pedigree_info{mat_id, am_dist, pat_id, ap_dist, p_score}
						candidate_pedigrees_info[i] = p_info
						//	fmt.Println(" pedigrees info:  ", p_info) // candidate_pedigrees_info[i])
						i++
					}
				}

				sort.Slice(candidate_pedigrees_info,
					func(i, j int) bool {
						return candidate_pedigrees_info[i].Amp_pedigree_score < candidate_pedigrees_info[j].Amp_pedigree_score
					})

				if pas, ok := pedigrees_and_scores[acc_id]; ok {
					fmt.Printf("%s  %s %6.4f %s %6.4f %6.4f  ", acc_id, pas.Mat_id, pas.Am_dist, pas.Pat_id, pas.Ap_dist, pas.Amp_pedigree_score)
				} else {
					fmt.Printf("%s %s  ", acc_id, "X -1 Y -1 -1")
				}

				for j, pas := range candidate_pedigrees_info {
					fmt.Printf("%s %6.4f %s %6.4f %6.4f  ", pas.Mat_id, pas.Am_dist, pas.Pat_id, pas.Ap_dist, pas.Amp_pedigree_score)
					if j >= 5 {
						break
					}
				}
				fmt.Printf("\n")
				
			}
		//	line_number++
		}
	}

}

// ********************************************************************************

// returns g_count, b_count pedigree_score = b_count/(g_count + b_count)
func calculate_pedigree_score(a_gt string, m_gt string, p_gt string) (int, int) {
	nmarkers := len(a_gt)
	if len(m_gt) != nmarkers || len(p_gt) != nmarkers {
		os.Exit(1)
	}
	var g_count, b_count = 0, 0
	for i := 0; i < len(a_gt); i++ {
		//	c1 := seq1[i : i+1]
		a := a_gt[i : i+1]
		m := m_gt[i : i+1]
		p := p_gt[i : i+1]

		// var g0022, b0022, g0110, b0110, g1221, b1221, g0220, b0220, g11 = 0, 0, 0, 0, 0, 0, 0, 0, 0
		if m == "0" {
			if p == "0" {
				if a == "0" {
					g_count++
				} else {
					b_count++
				}
			} else if p == "1" {
				if a == "2" {
					b_count++
				} else {
					g_count++
				}
			} else { // p == 2
				if a == "1" {
					g_count++
				} else {
					b_count++
				}
			}
		} else if m == "1" {
			if p == "0" {
				if a == "2" {
					b_count++
				} else {
					g_count++
				}
			} else if p == "1" {
				//	g_count++ // or don't increment
			} else { // p == 2
				if a == "0" {
					b_count++
				} else {
					g_count++
				}
			}
		} else if m == "2" {
			if p == "0" {
				if a == "1" {
					g_count++
				} else {
					b_count++
				}
			} else if p == "1" {
				if a == "0" {
					b_count++
				} else {
					g_count++
				}
			} else { // p == 2
				if a == "2" {
					g_count++
				} else {
					b_count++
				}
			}
		}
	}
	return g_count, b_count
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

func calculate_candidate_distances(id_seqset map[string]*sequenceset.Sequence_set, qid_cmfpq map[string]*priorityqueue.PriorityQueue) int {
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
			dist, ok := idpair_dist[idpair] /* */
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

func Initialize_priorityqueues(seq_set *sequenceset.Sequence_set, qid_cmfpq *map[string]*priorityqueue.PriorityQueue) {
	// make an empty priority queue for each sequence to hold the best candidate matches to that sequence
	for qid, _ := range seq_set.SeqId_index {
		pq := make(priorityqueue.PriorityQueue, 0)
		heap.Init(&pq)
		(*qid_cmfpq)[qid] = &pq
	}
}

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

func store_matches(ch chan map[string][]*mytypes.MatchInfo, qid_allokmatches map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
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

}

// get the top n_keep matches to each q seq
func search(qscs *seqchunkset.Sequence_chunk_set, scs *seqchunkset.Sequence_chunk_set, n_keep int, ch chan map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
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
}

func send_top_candidates(qid_matchinfos map[string][]*mytypes.MatchInfo, id_seqset map[string]*sequenceset.Sequence_set, n_keep int, wg *sync.WaitGroup,
	ch chan []*mytypes.IdSeq) {
	defer wg.Done()
	for q_id, okmatches := range qid_matchinfos { // for each query there are n_cpus*n_keep candidates
		sort.Slice(okmatches, func(i, j int) bool { // sort them by ChunkMatchFraction
			return okmatches[i].ChunkMatchFraction > okmatches[j].ChunkMatchFraction
		}) /* */

		qss := id_seqset[q_id]
		//	_ = qss.Check_seq_index_id_maps()
		q_seq := qss.Sequences[qss.SeqId_index[q_id]]
		//	top_matches := make([]mytypes.IdCmfDistance, n_keep)

		q_idseq := mytypes.IdSeq{q_id, q_seq}

		n_to_do := etc.MinInt(n_keep, len(okmatches))
		/*	if len(okmatches) < n_keep {
			n_to_do = len(okmatches)
		} /* */

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
	/* for {
	q_and_matches, ok := <-ch
	if !ok {
		break
	} */
	for q_and_matches := range ch {
		query := q_and_matches[0]
		fmt.Printf("%s  ", query.Id)
		for i := 1; i < len(q_and_matches); i++ {
			subj := q_and_matches[i]
			fmt.Printf("%s %7.5f  ", subj.Id, subj.Distance)
		}
		fmt.Println("")
	}
}

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
