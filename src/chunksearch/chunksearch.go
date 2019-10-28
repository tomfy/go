package main

import (
	//	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	//	"regexp"
	//	"strings"
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
	chunk_size := flag.Int("size", 6, "number of snps in each chunk")
	n_chunks := flag.Int("chunks", -1, "number of chunks to use")
	n_keep := flag.Int("keep", 20, "# of best matches to keep")

	flag.Parse()

	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	fmt.Printf("# file1: %s\n", *file1)
	fmt.Printf("# file2: %s\n", *file2)
	fmt.Printf("# chunk size: %d\n", *chunk_size)
	fmt.Printf("# keep: %d\n", *n_keep)

	t0 := time.Now()
	sequence_set1 := sequenceset.Construct_from_fasta_file(*file1)
	fmt.Fprintf(os.Stderr, "Done constructing sequence set 1.\n")
fmt.Fprintf(os.Stderr, "chunk size: %d  n_chunks: %d \n", *chunk_size, *n_chunks)
	seqchset := seqchunkset.Construct_from_sequence_set(sequence_set1, *chunk_size, *n_chunks)
	fmt.Fprintf(os.Stderr, "Done constructing sequence chunk set.\n")
	n_seqs1 := len(seqchset.Sequence_set.Sequences)
	t1 := time.Now()
	fmt.Fprintf(os.Stderr, "# time to construct: %v \n", t1.Sub(t0))
	fmt.Printf("# time to construct: %v \n", t1.Sub(t0))

	sequence_set2 := sequenceset.Construct_from_fasta_file(*file2)
	fmt.Fprintf(os.Stderr, "Done constructing sequence set 2.\n")

	// for each sequence in set 2, search for relatives in set 1

	for index, seq2 := range sequence_set2.Sequences {
		//	s2 := sequence_set2.Seqs[index]
		id2 := sequence_set2.Index_to_id(index)
		mindex_count_pairs := make([]mytypes.Pair_int_int, n_seqs1)
		for i := 0; i < n_seqs1; i++ {
			mindex_count_pairs[i].A = i
		}

		top_mindex_count_pairs := seqchset.Get_chunk_matchindex_counts(seq2, mindex_count_pairs, *n_keep)
		if index%500 == 0 {
			fmt.Fprintf(os.Stderr, "Search %d done.\n", index)
		}
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
	t2 := time.Now()

	fmt.Fprintf(os.Stderr, "All searches done.\n")
	fmt.Printf("# time to search: %v \n", t2.Sub(t1))
}

func distance_x(seq1 []uint, seq2 []uint) float64 {
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
