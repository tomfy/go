package main

import (
	//	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	//	"regexp"
	//	"strings"
	"math/rand"
	"mytypes"
	"runtime/pprof"
	"seqchunkset"
	"sequenceset"
	"sort"
	"time"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpy profile to file")

// read in one set of sequences, and do all N choose 2 chunkwise comparisons,
// read a sequence in, compare to all the previous ones, store it.

func main() {

	/* command line options: */

	/* input file: */
	var file string
	flag.StringVar(&file, "f1", "", "name of first fasta file.")

	/* search control parameters */
	var chunk_size, n_chunks, n_keep, n_reps  int
	var seed int64
	flag.IntVar(&chunk_size, "size", 5, "number of snps in each chunk")
	flag.IntVar(&n_chunks, "chunks", -1, "number of chunks to use")
	flag.IntVar(&n_keep, "keep", 20, "# of best matches to keep")
	flag.Int64Var(&seed, "seed", -1, "# rng seed (default: set from clock.)")
	flag.IntVar(&n_reps, "reps", 1, "# number of times to repeat the whole search (with different random chunk sets)")
	var missing_data_prob float64
	flag.Float64Var(&missing_data_prob, "miss", -1, "# fraction missing data in genotypes")

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

	for irep := 0; irep < n_reps; irep++ {
		t0 := time.Now()
		sequence_set := sequenceset.Construct_from_fasta_file(file)
		sequence_set.Add_missing_data(missing_data_prob)
		t1 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to construct sequence set: %v \n", t1.Sub(t0))

		if n_chunks < 0 {
			n_chunks = int(sequence_set.Sequence_length / chunk_size)
		}
		fmt.Fprintf(os.Stderr, "# file: %s   ", file)
		fmt.Fprintf(os.Stderr, "# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)
		fmt.Printf("# file: %s    ", file)
		fmt.Printf("# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)

		// for each sequence in set:
		// search for related sequences among those that have been stored previously;
		// then store latest sequence.
		t2 := time.Now()
		seqchset := seqchunkset.Construct_empty(sequence_set.Sequence_length, chunk_size, n_chunks) // 
		id2__index1_matchcount := make(map[string][]*mytypes.IntIntIntF64) // keys strings (id2), values: slices
		for index2, seq2 := range sequence_set.Sequences {
			id2 := sequence_set.Index_to_id(index2)
			if index2 > 0 { // search against the previously read in sequences
				top_chunkwise_matches := seqchset.Get_chunk_matchindex_counts(seq2, n_keep)
				if index2%1000 == 0 {
					fmt.Fprintf(os.Stderr, "Search %d done.\n", index2)
				}
				id2__index1_matchcount[id2] = top_chunkwise_matches
			}			
			seqchset.Add_sequence(id2, seq2) // add latest sequence 
		}
		t3 := time.Now()
		fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
		fmt.Fprintf(os.Stderr, "# time to search: %v \n", t3.Sub(t2))

		for index2, seq2 := range sequence_set.Sequences {
			id2 := sequence_set.Index_to_id(index2)
			top_mindex_count_pairs := id2__index1_matchcount[id2]
			
			id_matchcount_distance_triples := make([]mytypes.Triple_string_int_double, len(top_mindex_count_pairs) )
			fmt.Printf("%s   ", id2)
			for i, mcp := range top_mindex_count_pairs {
				seq1_index := mcp.A
				seq1_id := sequence_set.Index_to_id(seq1_index)
				seq1 := sequence_set.Sequences[seq1_index]
				dist12 := distance(seq1, seq2)
				//	fmt.Printf("# dists: %10.6f  %10.6f \n", dist12, disty12)
				id_matchcount_distance_triples[i] = mytypes.Triple_string_int_double{seq1_id, mcp.B, dist12}
			}
			sort.Slice(id_matchcount_distance_triples,
				func(i, j int) bool { return id_matchcount_distance_triples[i].C < id_matchcount_distance_triples[j].C })

			for _, a_triple := range id_matchcount_distance_triples {
				fmt.Printf("%s %d %6.5f  ", a_triple.A, a_triple.B, a_triple.C)
			}
			fmt.Printf("\n")
		}
		t4 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to calculate distances: %v \n", t4.Sub(t3))

		fmt.Printf("# time to construct sequence set: %v \n", t1.Sub(t0))
		fmt.Printf("# time to search for candidates: %v \n", t3.Sub(t2))
		fmt.Printf("# time to calculate distances: %v \n", t4.Sub(t3))
		fmt.Printf("# total time: %v \n", t4.Sub(t0))
	}
	tend := time.Now()
	fmt.Printf("# total time for %d reps: %v\n", n_reps, tend.Sub(tstart))
}

func distance(seq1 string, seq2 string) float64 {
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
}


/* func distance_y(seq1 string, seq2 string) float64 {
	ok_count := 0   // counts sites where neither seq has missing data
	dist_count := 0 // sums differences, i.e. 0-1 -> +=1, 0-2 -> += 2, ...
	for i := 0; i < len(seq1); i++ {
		c1 := seq1[i : i+1]
		c2 := seq2[i : i+1]
		if c1 == "0" {
			if c2 == "0" {
				ok_count++
			} else if c2 == "1" {
				ok_count++
				dist_count++
			} else if c2 == "2" {
				ok_count++
				dist_count += 2
			}
		} else if c1 == "1" {
			if c2 == "0" {
				ok_count++
				dist_count++
			} else if c2 == "1" {
				ok_count++
			} else if c2 == "2" {
				ok_count++
				dist_count++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				ok_count++
				dist_count += 2
			} else if c2 == "1" {
				ok_count++
				dist_count++
			} else if c2 == "2" {
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
} */

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
