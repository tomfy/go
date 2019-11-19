package main

import (
	"flag"
	"fmt"
	"log"
	"math/rand"
	"mytypes"
	"os"
	"runtime"
	"runtime/pprof"
	"seqchunkset"
	"sequenceset"
	"strings"
	//	"sort"
	"time"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpy profile to file")

// read in one set of sequences, and do all N choose 2 chunkwise comparisons,
// read a sequence in, compare to all the previous ones, store it.

func main() {

	/* command line options: */

	/* input file: */
	var files_string string
	flag.StringVar(&files_string, "f", "", "comma separated names of input fasta files.")

	/* search control parameters */
	var chunk_size, n_chunks, n_keep, n_reps int
	var seed int64
	flag.IntVar(&chunk_size, "size", 4, "number of snps in each chunk")
	flag.IntVar(&n_chunks, "chunks", -1, "number of chunks to use")
	flag.IntVar(&n_keep, "keep", 20, "# of best matches to keep")
	flag.Int64Var(&seed, "seed", -1, "# rng seed (default: set from clock.)")
	flag.IntVar(&n_reps, "reps", 1, "# number of times to repeat the whole search (with different random chunk sets)")
	var missing_data_prob, max_missing_data_proportion float64
	flag.Float64Var(&missing_data_prob, "miss", -1, "# fraction missing data in genotypes")
	flag.Float64Var(&max_missing_data_proportion, "max_md", 0.1, "# max proportion of missing data to use snp in chunk set")

	flag.Parse()

	tstart := time.Now()
	if seed < 0 {
		seed = int64(tstart.Nanosecond())
	}
	rand.Seed(seed)

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}
	//

	files := strings.Split(files_string, ",")

//	search_time := time.Now()
//	distance_time := search_time
	search_time := int64(0)
	distance_time := int64(0)
	for irep := 0; irep < n_reps; irep++ {
		data_sets := make([]*seqchunkset.Sequence_chunk_set, 0, len(files))
		cumulative_total_chunk_match_count := 0
		cumulative_total_mdmd_match_count := 0
		qid_matches := make(map[string][]mytypes.StringF64F64) // keys: query ids, values: slices of {subj_id, chunkmatchfraction , dist}
		for _, file := range files {
		
			q_sequence_set := sequenceset.Construct_from_fasta_file(file, max_missing_data_proportion, missing_data_prob)
			// sequence_set.Add_missing_data(missing_data_prob)

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


			//
			fmt.Println("n data sets: ", len(data_sets))
			for _, data_set := range data_sets{
				t_before := time.Now()
				qid_matchcandidates, total_chunk_match_count, total_mdmd_match_count := data_set.Search(q_sequence_set, n_keep)
			//	fmt.Println("# XXXX search time: ", time.Now().Sub(t_before))
				search_time += int64(time.Now().Sub(t_before))
				cumulative_total_chunk_match_count += total_chunk_match_count
				cumulative_total_mdmd_match_count += total_mdmd_match_count
				//	fmt.Fprintln(os.Stdout, len(qid_matchcandidates), total_chunk_match_count, total_mdmd_match_count)
				t_before = time.Now()
				q_sequence_set.Candidate_distances_AB(data_set.Sequence_set, qid_matchcandidates, qid_matches)
				distance_time += int64(time.Now().Sub(t_before))
			}

			t_before := time.Now()
			seqchset := seqchunkset.Construct_empty(q_sequence_set, chunk_size, n_chunks) //
			qid_matchcandidates, total_chunk_match_count, total_mdmd_match_count := seqchset.Search_and_construct(n_keep)
			data_sets = append(data_sets, seqchset)
			search_time += int64(time.Now().Sub(t_before))
			fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
			fmt.Fprintf(os.Stderr, "# chunk match counts; neither md: %d, both md: %d\n",
				total_chunk_match_count, total_mdmd_match_count)
			fmt.Fprintln(os.Stderr, MemUsageString())
  
			//	q_sequence_set.Candidate_distances_AA(qid_matchcandidates, qid_matches)
			t_before = time.Now()
			q_sequence_set.Candidate_distances_AB(q_sequence_set, qid_matchcandidates, qid_matches)
			distance_time += int64(time.Now().Sub(t_before))
			memstring := MemUsageString()
			fmt.Fprintln(os.Stderr, memstring)
			fmt.Printf("# chunk match counts; neither md: %d, both md: %d\n", total_chunk_match_count, total_mdmd_match_count)
			fmt.Println(memstring)

		} // end loop over input files (data sets)

	} // end loop over reps

	tend := time.Now()
	fmt.Printf("# total time for %d reps: %v\n", n_reps, tend.Sub(tstart))
	fmt.Printf("# search time: %10.3f  distance_time: %10.3f \n", 0.001*float64(search_time/1000000), 0.001*float64(distance_time/1000000))
}

// ******************************************************************************

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


