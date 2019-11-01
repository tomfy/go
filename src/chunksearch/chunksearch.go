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

	file1 := flag.String("f1", "", "name of first fasta file.")
	file2 := flag.String("f2", "", "name of 2nd fasta file.")

	/* search control parameters */
	chunk_size := flag.Int("size", 5, "number of snps in each chunk")
	n_chunks := flag.Int("chunks", -1, "number of chunks to use")
	n_keep := flag.Int("keep", 20, "# of best matches to keep")
	seed_ptr := flag.Int("seed", -1, "# rng seed (default: set from clock.)")
	n_reps := flag.Int("reps", 1, "# number of times to repeat the whole search (with different random chunk sets)")
	missing_data_prob := flag.Float64("miss", -1, "# fraction missing data in genotypes")

	seed := int64(*seed_ptr)
	tstart := time.Now()
	if seed < 0 {
		seed = int64(tstart.Nanosecond())
	}
	rand.Seed(seed)

	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	fmt.Printf("# file1: %s   file2: %s\n", *file1, *file2)
	fmt.Printf("# chunk_size: %d   n_chunks: %d   n_keep: %d   seed: %d\n", *chunk_size, *n_chunks, *n_keep, seed)

	for irep := 0; irep < *n_reps; irep++ {
		t0 := time.Now()
		sequence_set1 := sequenceset.Construct_from_fasta_file(*file1)
		sequence_set1.Add_missing_data(*missing_data_prob)
		fmt.Fprintf(os.Stderr, "Done constructing sequence set 1.\n")
		fmt.Fprintf(os.Stderr, "chunk size: %d  n_chunks: %d  n_keep: %d\n", *chunk_size, *n_chunks, *n_keep)
		seqchset := seqchunkset.Construct_from_sequence_set(sequence_set1, *chunk_size, *n_chunks)
		fmt.Fprintf(os.Stderr, "Done constructing sequence chunk set.\n")
		n_seqs1 := len(seqchset.Sequence_set.Sequences)
		t1 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to construct: %v \n", t1.Sub(t0))

		sequence_set2 := sequenceset.Construct_from_fasta_file(*file2)
		sequence_set2.Add_missing_data(*missing_data_prob)
		fmt.Fprintf(os.Stderr, "# Done constructing sequence set 2.\n")

		// for each sequence in set 2, search for relatives in set 1

		id2__index1_matchcount := make(map[string][]mytypes.Pair_int_int) // keys strings (id2), values: slices
		for index2, seq2 := range sequence_set2.Sequences {
			//	s2 := sequence_set2.Seqs[index]
			id2 := sequence_set2.Index_to_id(index2)
			mindex_count_pairs := make([]mytypes.Pair_int_int, n_seqs1)
			for i := 0; i < n_seqs1; i++ {
				mindex_count_pairs[i].A = i
			}

			top_mindex_count_pairs := seqchset.Get_chunk_matchindex_counts(seq2, mindex_count_pairs, *n_keep)
			if index2%1000 == 0 {
				fmt.Fprintf(os.Stderr, "Search %d done.\n", index2)
			}
			id2__index1_matchcount[id2] = top_mindex_count_pairs
		}
		t2 := time.Now()
		fmt.Fprintf(os.Stderr, "# All searches done.\n")
		fmt.Fprintf(os.Stderr, "# time to search: %v \n", t2.Sub(t1))

		for index2, seq2 := range sequence_set2.Sequences {
			//	for id2, top_mindex_count_pairs := range id2__index1_matchcount {
			id2 := sequence_set2.Index_to_id(index2)
			top_mindex_count_pairs := id2__index1_matchcount[id2]

			id_matchcount_distance_triples := make([]mytypes.Triple_string_int_double, *n_keep, *n_keep)
			fmt.Printf("%s   ", id2)
			for i, mcp := range top_mindex_count_pairs {
				seq1_index := mcp.A
				seq1_id := sequence_set1.Index_to_id(seq1_index)
				seq1 := sequence_set1.Sequences[seq1_index]
				dist12 := distance(seq1, seq2)
				//	s1 := sequence_set1.Seqs[seq1_index]
				//	dist12_x := distance_x(s1.Gs, s2.Gs)
				id_matchcount_distance_triples[i] = mytypes.Triple_string_int_double{seq1_id, mcp.B, dist12}
				//	fmt.Printf("%s %d %6.5f  ", seq1_id, mcp.B, dist12)
			}
			sort.Slice(id_matchcount_distance_triples, func(i, j int) bool { return id_matchcount_distance_triples[i].C < id_matchcount_distance_triples[j].C })

			for _, a_triple := range id_matchcount_distance_triples {
				fmt.Printf("%s %d %6.5f  ", a_triple.A, a_triple.B, a_triple.C)
			}
			fmt.Printf("\n")
		}
		t3 := time.Now()
		fmt.Fprintf(os.Stderr, "# time to calculate distances: %v \n", t3.Sub(t2))

		fmt.Printf("# time to construct: %v \n", t1.Sub(t0))
		fmt.Printf("# time to search: %v \n", t2.Sub(t1))
		fmt.Printf("# time to calculate distances: %v \n", t3.Sub(t2))
		fmt.Printf("# total time: %v \n", t3.Sub(t0))
	}
	tend := time.Now()
	fmt.Printf("# total time for %d reps: %v\n", *n_reps, tend.Sub(tstart))
}

func distance(seq1 string, seq2 string) float64 {
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
