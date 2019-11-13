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
	flag.StringVar(&file, "f", "", "name of first fasta file.")

	/* search control parameters */
	var chunk_size, n_chunks, n_keep, n_reps int
	var seed int64
	flag.IntVar(&chunk_size, "size", 4, "number of snps in each chunk")
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
		sequence_set := sequenceset.Construct_from_fasta_file(file, missing_data_prob)
		// sequence_set.Add_missing_data(missing_data_prob)
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
		// search for candidate related sequences among those that have been stored previously;
		// then store latest sequence.
		t2 := time.Now()
		seqchset := seqchunkset.Construct_empty(sequence_set.Sequence_length, chunk_size, n_chunks) //
		qid_smatchinfos := make(map[string][]*mytypes.IntIntIntF64)                                 // keys strings (id2), values: slices
		for qindex, qseq := range sequence_set.Sequences {
			qid := sequence_set.Seq_index_to_id(qindex)
			if qindex > 0 { // search against the previously read in sequences
				top_smatchinfos := seqchset.Get_chunk_matchindex_counts(qseq, n_keep)
				if qindex%1000 == 0 {
					fmt.Fprintf(os.Stderr, "Search %d done.\n", qindex)
				}
				qid_smatchinfos[qid] = top_smatchinfos
			}
			seqchset.Add_sequence(qid, qseq) // add latest sequence
		}
		t3 := time.Now()
		fmt.Fprintf(os.Stderr, "# All searches for candidates done.\n")
		fmt.Fprintf(os.Stderr, "# time to search: %v \n", t3.Sub(t2))

		fmt.Fprintln(os.Stderr, MemUsageString())
		// do the full distance calculation for each candidate found above based on chunkwise analysis
		for qindex, qseq := range sequence_set.Sequences {
			qid := sequence_set.Seq_index_to_id(qindex)
			top_smatchinfos := qid_smatchinfos[qid]

			id_matchcount_distance_triples := make([]mytypes.Triple_string_int_double, len(top_smatchinfos))
			fmt.Printf("%s   ", qid)
			for i, smatchinfo := range top_smatchinfos {
				sseq_index := smatchinfo.A
				sseq_id := sequence_set.Seq_index_to_id(sseq_index)
				sseq := sequence_set.Sequences[sseq_index]
			//	dist_old := distance_old(sseq, qseq)
				n00_22, n11, nd1, nd2 := distance(sseq, qseq)
				dist := float64(nd1 + 2*nd2)/float64(n00_22 + n11 + nd1 + nd2)
			//	fmt.Printf("%v  %v\n", dist_old, dist)
				/*if(dist != distx){
					os.Exit(1)
				}*/
				id_matchcount_distance_triples[i] = mytypes.Triple_string_int_double{sseq_id, smatchinfo.B, dist}
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
		memstring := MemUsageString()
		fmt.Fprintln(os.Stderr, memstring)

		fmt.Printf("# time to construct sequence set: %v \n", t1.Sub(t0))
		fmt.Printf("# time to search for candidates: %v \n", t3.Sub(t2))
		fmt.Printf("# time to calculate distances: %v \n", t4.Sub(t3))
		fmt.Printf("# total time: %v \n", t4.Sub(t0))
		fmt.Println(memstring)
	}
	tend := time.Now()
	fmt.Printf("# total time for %d reps: %v\n", n_reps, tend.Sub(tstart))
}

func distance_old(seq1 string, seq2 string) float64 {
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

func distance(seq1 string, seq2 string) (int, int, int, int) {
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
	/*	ok_count := n00_22 + n11 + one_count + two_count // number of sites where neither seq has missing data
		dist_count := one_count + 2*two_count          // sums differences, i.e. 0-1 -> +=1, 0-2 -> += 2, ...
		var distance float64
		if ok_count > 0 {
			distance = float64(dist_count) / float64(ok_count)
		} else {
			distance = -1.0 // couldn't calculate because no sites without missing data
		}
		return distance */
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
