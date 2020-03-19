package main

import (
	"bufio"
	"container/heap"
	"etc"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strings"
	"sync"
	"time"

	"mytypes"
	"priorityqueue"
	"seqchunkset"
	"sequenceset"
)

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
	var n_keep, n_cpus int
	var pauses bool
	var seed int64
	max_md_prop := 0
	flag.IntVar(&n_cpus, "cpus", 1, "Number of cpus to use.")
	flag.IntVar(&n_keep, "keep", 20, "Number of best matches to keep.")
	flag.Int64Var(&seed, "seed", -1, "Rng seed (default: set from clock.)")
	/*	flag.IntVar(&n_reps, "reps", 1, "Number of times to repeat the whole search with different random chunk sets.")
		var missing_data_prob, max_missing_data_proportion float64
		flag.Float64Var(&missing_data_prob, "miss", -1, "Proportion missing data to add to genotypes.")
		flag.Float64Var(&max_missing_data_proportion, "max_md", 0.1, "Max proportion of missing data to use marker in chunk set.")
		flag.BoolVar(&pauses, "pause", false, "false: no pauses, true: pause at various points until enter key pressed.") /* */
	flag.Parse() // parse the command line

	// Seed the rng
	if seed < 0 {
		seed = int64(time.Now().Nanosecond())
	}
	rand.Seed(seed)

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

	var seq_set *sequenceset.Sequence_set
	sequenceset.Construct_sets_from_matrix_file(gt_matrix_file, n_cpus, max_missing_data_proportion, &id_seqset, seq_set)
	seq_set := Construct_from_matrix_file(gt_matrix_file, max_md_prop, id_seqset, seq_set, wg1)

	n_markers := len(seq_sets[0].Sequences[0])
	fmt.Fprintln(os.Stderr, "# n markers: ", n_markers)

	//****************** load closest relatives from file (if specified):

	//****************** read pedigree information from file:

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
			acc_gt := seq_set.Sequences[seq_set.SeqId_index[acc_id]]
			mat_gt := seq_set.Sequences[seq_set.SeqId_index[mat_id]]
			pat_gt := seq_set.Sequences[seq_set.SeqId_index[pat_id]]
			g_count, b_count := calculate_pedigree_score(acc_gt, mat_gt, pat_gr)
			pedigree_score := float64(b_count) / float64(g_count+b_count)
			fmt.Println(acc_id, mat_id, pat_id, pedigree_score)
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
	for i := 0; i < len(a_gt); i++ {
		a := a_gt[i]
		m := m_gt[i]
		p := p_gt[i]

		var g_count, b_count = 0, 0
		// var g0022, b0022, g0110, b0110, g1221, b1221, g0220, b0220, g11 = 0, 0, 0, 0, 0, 0, 0, 0, 0
		if m == 0 {
			if p == 0 {
				if a == 0 {
					g_count++
				} else {
					b_count++
				}
			} else if p == 1 {
				if a == 2 {
					b_count++
				} else {
					g_count++
				}
			} else { // p == 2
				if a == 1 {
					g_count++
				} else {
					b_count++
				}
			}
		} else if m == 1 {
			if p == 0 {
				if a == 2 {
					b_count++
				} else {
					g_count++
				}
			} else if p == 1 {
				//	g_count++ // or don't increment
			} else { // p == 2
				if a == 0 {
					b_count++
				} else {
					g_count++
				}
			}
		} else if m == 2 {
			if p == 0 {
				if a == 1 {
					g_count++
				} else {
					b_count++
				}
			} else if p == 1 {
				if a == 0 {
					b_count++
				} else {
					g_count++
				}
			} else { // p == 2
				if a == 2 {
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
func calculate_distances(inch chan []*mytypes.IdSeq, outch chan []*mytypes.IdCmfDistance, wg *sync.WaitGroup) {
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
			n00_22, n11, nd1, nd2 := sequenceset.Distance(query.Sequence, subj.Sequence)
			dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
			sout := mytypes.IdCmfDistance{subj.Id, -1, dist}
			out[i] = &sout
		}
		sort.Slice(out, func(i, j int) bool { // sort them by ChunkMatchFraction
			return out[i].Distance < out[j].Distance
		}) /* */
		out[0].Distance = 0 // now set it to 0
		outch <- out
	}
}

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
