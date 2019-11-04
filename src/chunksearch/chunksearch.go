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

func main() {

	/* command line options: */

	/* input files: */

	file1_ptr := flag.String("f1", "", "name of first fasta file.")
	file2_ptr := flag.String("f2", "", "name of 2nd fasta file.")

	/* search control parameters */
	chunk_size_ptr := flag.Int("size", 5, "number of snps in each chunk")
	n_chunks_ptr := flag.Int("chunks", -1, "number of chunks to use")
	n_keep_ptr := flag.Int("keep", 20, "# of best matches to keep")
	seed_ptr := flag.Int("seed", -1, "# rng seed (default: set from clock.)")
	n_reps_ptr := flag.Int("reps", 1, "# number of times to repeat the whole search (with different random chunk sets)")
	missing_data_prob_ptr := flag.Float64("miss", -1, "# fraction missing data in genotypes")

	flag.Parse()

	file1 := *file1_ptr
	file2 := *file2_ptr
	chunk_size := *chunk_size_ptr
	n_chunks := *n_chunks_ptr
	n_keep := *n_keep_ptr
	n_reps := *n_reps_ptr
	missing_data_prob := *missing_data_prob_ptr

	seed := int64(*seed_ptr)
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
		sequence_set1 := sequenceset.Construct_from_fasta_file(file1)
		sequence_set1.Add_missing_data(missing_data_prob)
		t1 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to construct sequence set 1: %v \n", t1.Sub(t0))

		if n_chunks < 0 {
			n_chunks = int(sequence_set1.Sequence_length / chunk_size)
		}
		fmt.Fprintf(os.Stderr, "# file1: %s   file2: %s\n", file1, file2)
		fmt.Fprintf(os.Stderr, "# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)
		fmt.Printf("# file1: %s   file2: %s\n", file1, file2)
		fmt.Printf("# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", chunk_size, n_chunks, n_keep, seed)

		seqchset := seqchunkset.Construct_from_sequence_set(sequence_set1, chunk_size, n_chunks)
		t2 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to construct sequence chunk set: %v \n", t2.Sub(t1))

		sequence_set2 := sequenceset.Construct_from_fasta_file(file2)
		sequence_set2.Add_missing_data(missing_data_prob)
		t3 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to construct sequence set 2: %v \n", t3.Sub(t2))

		// for each sequence in set 2, search for relatives in set 1

		id2__index1_matchcount := make(map[string][]*mytypes.IntIntIntF64) // keys strings (id2), values: slices
		for index2, seq2 := range sequence_set2.Sequences {
			id2 := sequence_set2.Index_to_id(index2)
			top_chunkwise_matches := seqchset.Get_chunk_matchindex_counts(seq2, n_keep)
			if index2%1000 == 0 {
				fmt.Fprintf(os.Stderr, "Search %d done.\n", index2)
			}
			id2__index1_matchcount[id2] = top_chunkwise_matches
		}
		t4 := time.Now()
		fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
		fmt.Fprintf(os.Stderr, "# time to search: %v \n", t4.Sub(t3))

		for index2, seq2 := range sequence_set2.Sequences {
			id2 := sequence_set2.Index_to_id(index2)
			top_mindex_count_pairs := id2__index1_matchcount[id2]

			id_matchcount_distance_triples := make([]mytypes.Triple_string_int_double, n_keep)
			fmt.Printf("%s   ", id2)
			for i, mcp := range top_mindex_count_pairs {
				seq1_index := mcp.A
				seq1_id := sequence_set1.Index_to_id(seq1_index)
				seq1 := sequence_set1.Sequences[seq1_index]
				//	disty12 := distance_y(seq1, seq2)
				dist12 := distance(seq1, seq2)
				//	fmt.Printf("# dists: %10.6f  %10.6f \n", dist12, disty12)
				id_matchcount_distance_triples[i] = mytypes.Triple_string_int_double{seq1_id, mcp.B, dist12}
			}
			sort.Slice(id_matchcount_distance_triples, func(i, j int) bool { return id_matchcount_distance_triples[i].C < id_matchcount_distance_triples[j].C })

			for _, a_triple := range id_matchcount_distance_triples {
				fmt.Printf("%s %d %6.5f  ", a_triple.A, a_triple.B, a_triple.C)
			}
			fmt.Printf("\n")
		}
		t5 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to calculate distances: %v \n", t5.Sub(t4))

		fmt.Printf("# time to construct sequence set 1: %v \n", t1.Sub(t0))
		fmt.Printf("# time to construct sequence chunk set: %v \n", t2.Sub(t1))
		fmt.Printf("# time to construct sequence set 2: %v \n", t3.Sub(t2))
		fmt.Printf("# time to search for candidates: %v \n", t4.Sub(t3))
		fmt.Printf("# time to calculate distances: %v \n", t5.Sub(t4))
		fmt.Printf("# total time: %v \n", t5.Sub(t0))
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
				//	dist_count++
			} else if c2 == "2" {
				two_count++
				//	dist_count += 2
			}
		} else if c1 == "1" {
			if c2 == "0" {
				one_count++
				//	dist_count++
			} else if c2 == "1" {
				zero_count++
			} else if c2 == "2" {
				one_count++
				//	dist_count++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				two_count++
				//	dist_count += 2
			} else if c2 == "1" {
				one_count++
				//	dist_count++
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

func distance_y(seq1 string, seq2 string) float64 {
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
}

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
