package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math/rand"
	"os"
	"runtime"
	"runtime/pprof"
	"sort"
	"strings"
	"time"
	//	"reflect"
	"container/heap"
	"sync"

	"etc"
	"mytypes"
	"priorityqueue"
	"seqchunkset"
	"sequenceset"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpy profile to file")

// read in one set of sequences, and do all N choose 2 chunkwise comparisons,
// read a sequence in, compare to all the previous ones, store it.

func main() {

	reader := bufio.NewReader(os.Stdin)
	/* command line options: */

	/* input file: */
	var files_string string
	flag.StringVar(&files_string, "f", "", "comma separated names of input fasta files.")

	/* search control parameters */
	var chunk_size, n_chunks, n_keep, n_reps, n_cpus, way int
	var seed int64
	flag.IntVar(&n_cpus, "cpus", 1, "n cpus (number of cpus to use)")
	flag.IntVar(&chunk_size, "size", 4, "number of snps in each chunk")
	flag.IntVar(&n_chunks, "chunks", -1, "number of chunks to use") // default of -1  ->  each marker used once.
	flag.IntVar(&way, "way", 1, "1: old way, other: new way.")
	flag.IntVar(&n_keep, "keep", 20, "# of best matches to keep")
	//	flag.IntVar(&way, "way", 1, "# 1 (default) old way, other: new way")
	flag.Int64Var(&seed, "seed", -1, "# rng seed (default: set from clock.)")
	flag.IntVar(&n_reps, "reps", 1, "# number of times to repeat the whole search (with different random chunk sets)")
	var missing_data_prob, max_missing_data_proportion float64
	flag.Float64Var(&missing_data_prob, "miss", -1, "# fraction missing data in genotypes")
	flag.Float64Var(&max_missing_data_proportion, "max_md", 0.1, "# max proportion of missing data to use snp in chunk set")
	//	flag.IntVar(&save2, "save2", 0, "# if 1, save factor or 2 by only doing 1 direction.")

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

	q_and_sfiles := strings.Split(files_string, ";") // split on semicolon

	setup_time := int64(0)
	q_setup_time := int64(0)
	s_setup_time := int64(0)
	search_time := int64(0)
	distance_time := int64(0)
	distance_calc_count := 0
	t_start := time.Now()

	//	n_cpus := 2

	var wg1 sync.WaitGroup
	var wg2 sync.WaitGroup
	for irep := 0; irep < n_reps; irep++ {
		t_start = time.Now()
		// can either do 1st file vs all others (just get best matches to sequences in first file) (mode 1)
		// or do all v all, getting best matches to every sequence (mode 2)
		if len(q_and_sfiles) > 1 { // 'mode 1'
			qfiles := strings.Split(q_and_sfiles[0], ",") // split on ,
			sfiles := strings.Split(q_and_sfiles[1], ",") // split on ,
			if len(qfiles) > 1 || len(sfiles) > 1 {
				fmt.Fprintln(os.Stderr, "only 1 qfile and 1 sfile allowed; will analyze 1st qfile v 1st sfile.")
				os.Exit(1)
			}
			qfile := qfiles[0]
			sfile := sfiles[0]

			id_seqset := make(map[string]*sequenceset.Sequence_set)

			//**************  setup query Sequence_set and Sequence_chunk_set :
			q_setup_start := time.Now()
			q_seq_sets := make([]*sequenceset.Sequence_set, n_cpus)
			
fmt.Fprintln(os.Stderr, "Before constructing q sequence sets. Press 'enter' to continue.")
			x, _ := reader.ReadString('\n')
			
			sequenceset.Construct_sets_from_matrix_file(qfile, n_cpus, max_missing_data_proportion, &id_seqset, q_seq_sets)

		
			fmt.Fprintln(os.Stderr, "After constructing q sequence sets. Press 'enter' to continue.")
			x, _ = reader.ReadString('\n')
			_ = x
			
			if n_chunks < 0 {
				n_chunks = int(q_seq_sets[0].Sequence_length / chunk_size)
			}
			n_markers := len(q_seq_sets[0].Sequences[0])
			fmt.Fprintln(os.Stderr, "# n markers: ", n_markers)
			the_chunk_specs := seqchunkset.Get_chunk_specs(q_seq_sets[0], chunk_size, n_chunks)
			fmt.Fprintln(os.Stderr, "After constructing chunk-specs. Press 'enter' to continue." )
			x, _ = reader.ReadString('\n')
			

			q_seqchunksets := seqchunkset.Construct_multiple_from_sequence_sets(q_seq_sets, the_chunk_specs, false, true) // chunk_size, n_chunks)
			q_setup_time = int64(time.Now().Sub(q_setup_start))

			
		
			fmt.Fprintln(os.Stderr, "After constructing q seqchunkset. Press 'enter' to continue: ")
			x, _ = reader.ReadString('\n')

			//***************  setup subj. Sequence_set and Sequence_chunk_set :
			s_setup_start := time.Now()
			s_seq_sets := make([]*sequenceset.Sequence_set, n_cpus)
			sequenceset.Construct_sets_from_matrix_file(sfile, n_cpus, max_missing_data_proportion, &id_seqset, s_seq_sets)
			//	fmt.Fprintln(os.Stderr, "Subj. len(s_seq_sets): ", len(s_seq_sets))
		
			fmt.Fprintln(os.Stderr, "After constructing subj. sequence set. Press 'enter' to continue.")
			x, _ = reader.ReadString('\n')
			
			s_seqchunksets := seqchunkset.Construct_multiple_from_sequence_sets(s_seq_sets, the_chunk_specs, false, true) // chunk_size, n_chunks)
			s_setup_time = int64(time.Now().Sub(s_setup_start))
			fmt.Fprintln(os.Stderr, "After constructing subj. seqchunkset. Press 'enter' to continue.")
			x, _ = reader.ReadString('\n')
			
			fmt.Fprintf(os.Stderr, "# chunk_size: %d  n_chunks: %d  n_keep: %d  n_cpus: %d  seed: %d\n",
				chunk_size, n_chunks, n_keep, n_cpus, seed)
			fmt.Fprintf(os.Stdout, "# chunk_size: %d  n_chunks: %d  n_keep: %d  n_cpus: %d  seed: %d\n",
				chunk_size, n_chunks, n_keep, n_cpus, seed)

			setup_time += int64(time.Now().Sub(t_start))
			//***************  end of setup

			cumulative_total_chunk_match_count := 0
			cumulative_total_mdmd_match_count := 0

			//****************  do the chunkwise search for candidates
			var qid_allokmatches map[string][]*mytypes.MatchInfo
			t_before := time.Now()
			qid_allokmatches = make(map[string][]*mytypes.MatchInfo)
			for j := 0; j < n_cpus; j++ { // j is the offset between q set index (i) and s set index (k)

				the_channel := make(chan map[string][]*mytypes.MatchInfo)
				wg1.Add(1)
				go store_matches(the_channel, qid_allokmatches, &wg1)

				for i, qscs := range q_seqchunksets {
					k := (i + j) % n_cpus
					wg2.Add(1)
					go search(qscs, s_seqchunksets[k], n_keep, the_channel, &wg2)
				}
				wg2.Wait()
				fmt.Fprintln(os.Stderr, "# Search with offset j: ", j, " done.")
				close(the_channel)
				wg1.Wait()
				fmt.Fprintln(os.Stderr, "After search with offset: ", j, ". Press 'enter' to continue.")
				x, _ = reader.ReadString('\n')
			}

			search_time += int64(time.Now().Sub(t_before))
			fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
			//****************  done with chunkwise (candidate) search

			fmt.Println("# n qids: ", len(qid_allokmatches))

			//****************  do the distance calculations for the top n_keep candidates
			t_before_dists := time.Now()
			channel_1 := make(chan []*mytypes.IdSeq)
			channel_2 := make(chan []*mytypes.IdCmfDistance)

			wg1.Add(1)
			go send_top_candidates(qid_allokmatches, id_seqset, n_keep, &wg1, channel_1)

			for i := 0; i < n_cpus; i++ {
				wg2.Add(1)
				go calculate_distances(channel_1, channel_2, &wg2)
			}
			fmt.Fprintln(os.Stderr, "Distance calculations underway. Press 'enter' to continue.")
				x, _ = reader.ReadString('\n')
			go output(channel_2)
			wg1.Wait()
			close(channel_1)
			wg2.Wait()
			close(channel_2)
			distance_time = int64(time.Now().Sub(t_before_dists))
				fmt.Fprintln(os.Stderr, "After distance calculations. Press 'enter' to continue.")
				x, _ = reader.ReadString('\n')
			//***************  distance calculations done

			fmt.Fprintf(os.Stderr, "# chunk match counts; neither md: %d, both md: %d\n",
				cumulative_total_chunk_match_count, cumulative_total_mdmd_match_count)
			fmt.Fprintln(os.Stderr, "# ", MemUsageString())

			memstring := MemUsageString()
			fmt.Println("# ", memstring)
			fmt.Printf("# chunk match counts; neither md: %d, both md: %d\n",
				cumulative_total_chunk_match_count, cumulative_total_mdmd_match_count)
			fmt.Println(memstring)

		} else { // MODE 2
			files := strings.Split(files_string, ",") // matrix files

			data_sets := make([]*seqchunkset.Sequence_chunk_set, 0, len(files))
			id_seqset := make(map[string]*sequenceset.Sequence_set)
			cumulative_total_chunk_match_count := 0
			cumulative_total_mdmd_match_count := 0

			qid_cmfpq := make(map[string]*priorityqueue.PriorityQueue) // one priority queue for each query
			for _, file := range files {

				q_sequence_set := &sequenceset.Sequence_set{}
				wg1.Add(1)
				go sequenceset.Construct_from_matrix_file(file, max_missing_data_proportion, &id_seqset, q_sequence_set, &wg1)
				wg1.Wait()

				Initialize_priorityqueues(q_sequence_set, &qid_cmfpq) // create an empty pq for each sequence in q_sequence_set, and add to qid_cmfpq

				if n_chunks < 0 {
					n_chunks = int(q_sequence_set.Sequence_length / chunk_size)
				}
				fmt.Fprintf(os.Stderr, "# file: %s   ", file)
				fmt.Fprintf(os.Stderr, "# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)
				fmt.Printf("# file: %s    ", file)
				fmt.Printf("# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)

				// for each sequence in set:
				// search for candidate related sequences among those that have been stored previously;
				// then store latest sequence

				fmt.Printf("# n data sets: %d\n", len(data_sets))
				for _, data_set := range data_sets {
					t_before := time.Now()
					qid_okmatches, qid_badmatches, total_chunk_match_count, total_mdmd_match_count := data_set.Search_pq(q_sequence_set, n_keep, false, &qid_cmfpq)
					_ = qid_okmatches

					if len(qid_badmatches) > 0 {
						fmt.Println("#  there are this many bad matches: ", len(qid_badmatches))
					}
					search_time += int64(time.Now().Sub(t_before))
					cumulative_total_chunk_match_count += total_chunk_match_count
					cumulative_total_mdmd_match_count += total_mdmd_match_count
				}

				t_before := time.Now()
				the_chunk_specs := seqchunkset.Get_chunk_specs(q_sequence_set, chunk_size, n_chunks)
				seqchset := seqchunkset.Construct_empty(q_sequence_set, the_chunk_specs) //
				//	fmt.Fprintln(os.Stderr, "save2: ", save2, (save2 == 1) )
				qid_okmatches, qid_badmatches, total_chunk_match_count, total_mdmd_match_count := seqchset.Search_pq(q_sequence_set, n_keep, true, &qid_cmfpq)
				_ = qid_okmatches
				if len(qid_badmatches) > 0 {
					fmt.Println("# there are bad matches.", len(qid_badmatches))
				}
				search_time += int64(time.Now().Sub(t_before))

				//	distance_calc_count += q_sequence_set.Candidate_distances_pq(q_sequence_set, qid_okmatches, qid_matches)

				data_sets = append(data_sets, seqchset)

				fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
				fmt.Fprintf(os.Stderr, "# chunk match counts; neither md: %d, both md: %d\n",
					total_chunk_match_count, total_mdmd_match_count)
				fmt.Fprintln(os.Stderr, "# ", MemUsageString())

				memstring := MemUsageString()
				fmt.Println("# ", memstring)
				fmt.Printf("# chunk match counts; neither md: %d, both md: %d\n", total_chunk_match_count, total_mdmd_match_count)
				fmt.Println(memstring)

			} // end loop over input files (data sets)
			t_before := time.Now()
			distance_calc_count += calculate_candidate_distances(id_seqset, qid_cmfpq) //, qid_matches)
			distance_time = int64(time.Now().Sub(t_before))                            /* */
		} // end mode 2
	} // end loop over reps

	tend := time.Now()
	fmt.Printf("# total time for %d reps: %v\n", n_reps, tend.Sub(t_start))
	fmt.Printf("# number of distances calculated: %d\n", distance_calc_count)
	fmt.Printf("# setup time(total, q, s): %6.3f %6.3f %6.3f  search time: %6.3f  distance_time: %6.3f \n",
		0.001*float64(setup_time/1000000), 0.001*float64(q_setup_time/1000000), 0.001*float64(s_setup_time/1000000),
		0.001*float64(search_time/1000000), 0.001*float64(distance_time/1000000))
}

//

// ******************************************************************************

// /*
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
			fmt.Printf("%s %5.3f %7.5f  ", x.Id, x.ChunkMatchFraction, x.Distance)
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

func store_matches(ch chan map[string][]*mytypes.MatchInfo, qid_allokmatches map[string][]*mytypes.MatchInfo, wg *sync.WaitGroup) {
	defer wg.Done()
	for {
		qid_okmatches, ok := <-ch // ok will be falst iff channel is empty and closed.
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
