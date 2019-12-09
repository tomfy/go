package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"math/rand"
	//	"mytypes"
	"os"
	//	"runtime"
	"runtime/pprof"
	"sort"
	"strings"
	"time"
	//
	//	"seqchunkset"
	"sequenceset"
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
	var dist_fraction float64
	var seed int64
	var top int
	flag.Int64Var(&seed, "seed", -1, "# rng seed (default: set from clock.)")
	flag.Float64Var(&dist_fraction, "dfrac", 1, "# fraction of all N choose 2 distances to do.")
	flag.IntVar(&top, "top", -1, "# output this many closest relatives to each.")

	max_missing_data_proportion := 0.1
	missing_data_prob := 0.0

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

	//	files := strings.Split(files_string, ",") // fasta files, or matrix files!

	/*	test_file := files[0] // filename of first file, now check to see whether is fasta or matrix
		input_format := what_is_format(test_file)
		fmt.Println("# input format: ", input_format)
		if input_format == "other" {
			os.Exit(1)
		} /* */

	var input_format string
	var qfiles, sfiles, files []string
	var qfile, sfile, file string
	var mode int
	var idpair_dist map[string]map[string]float64
	q_and_sfiles := strings.Split(files_string, ";") // split on ;
	if len(q_and_sfiles) > 1 {                       // 'mode 1'
		mode = 1
		qfiles = strings.Split(q_and_sfiles[0], ",") // split on ,
		sfiles = strings.Split(q_and_sfiles[1], ",") // split on ,
		input_format = what_is_format(sfiles[0])
		qfile = qfiles[0]
		sfile = sfiles[0]
		_ = qfile
		_ = sfile
	} else { // mode 2
		mode = 2
		files = strings.Split(files_string, ",") // fasta files, or matrix files!
		file = files[0]
		input_format = what_is_format(file)
		if len(files) > 1 {
			fmt.Fprintf(os.Stderr, "More than 1 file specified; will just use first one.\n")

		}

	}

	fmt.Fprintln(os.Stderr, "# mode: ", mode)
	fmt.Println("# input format: ", input_format)
	if input_format == "other" {
		os.Exit(1)
	}

	//	search_time := time.Now()
	//	distance_time := search_time

	//	for _, file := range files {

	var sequence_set1, sequence_set2 *sequenceset.Sequence_set
	if mode == 1 { // 1 query file and 1 subj. file
		id_seqset := make(map[string]*sequenceset.Sequence_set)
		if input_format == "fasta" {
			sequence_set1 = sequenceset.Construct_from_fasta_file(qfile, max_missing_data_proportion, missing_data_prob, &id_seqset)
			sequence_set2 = sequenceset.Construct_from_fasta_file(sfile, max_missing_data_proportion, missing_data_prob, &id_seqset)
		} else if input_format == "matrix" {
			sequence_set1 = sequenceset.Construct_from_matrix_file(qfile, max_missing_data_proportion, &id_seqset)
			sequence_set2 = sequenceset.Construct_from_matrix_file(sfile, max_missing_data_proportion, &id_seqset)
		}
		fmt.Fprintln(os.Stderr, len(sequence_set1.Sequences), len(sequence_set2.Sequences))
		idpair_dist = sequence_set1.Fraction_of_all_distances_AB(sequence_set2, dist_fraction)
	} else {
		id_seqset := make(map[string]*sequenceset.Sequence_set)
		if input_format == "fasta" {
			sequence_set1 = sequenceset.Construct_from_fasta_file(file, max_missing_data_proportion, missing_data_prob, &id_seqset)
		} else if input_format == "matrix" {
			sequence_set1 = sequenceset.Construct_from_matrix_file(file, max_missing_data_proportion, &id_seqset)
		}
		// sequence_set.Add_missing_data(missing_data_prob)
		idpair_dist = sequence_set1.Fraction_of_all_distances_AA(dist_fraction)
	}
	if top > 0 {
		fmt.Fprintln(os.Stderr, "n id1s: ", len(idpair_dist))
		for id1, id2_dist := range idpair_dist {
			sorted_keys := keys_sorted_by_value(id2_dist)
			fmt.Printf("%s  ", id1)
			for i2, id2 := range sorted_keys {
				d12 := id2_dist[id2]
				fmt.Printf("%s %8.5f  ", id2, d12)
				if i2 >= top-1 {
					break
				}
			}
			fmt.Println("")
		}
	} else {
		for id1, id2_dist := range idpair_dist {
			//	sorted_keys := keys_sorted_by_value(id2_dist)
			for id2, d12 := range id2_dist { // sorted_keys {
				// d12 := id2_dist[id2]
				fmt.Printf("%s %s %8.5f\n", id1, id2, d12)
			}
		}
	}
	// }
	tend := time.Now()
	fmt.Printf("time to execute: %v \n", tend.Sub(tstart))
}

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

// for a map w string keys, int values, return a slice containing the keys, but sorted with smallest value ones at beginning.
func keys_sorted_by_value(amap map[string]float64) []string {
	keys := make([]string, 0, len(amap))
	for k, _ := range amap {
		keys = append(keys, k)
		//	fmt.Fprintf(os.Stderr, "k v: %s  %d\n", k, v)
	}
	sort.Slice(keys, func(i, j int) bool { return amap[keys[i]] < amap[keys[j]] })
	//	fmt.Fprintf(os.Stderr, "max_mds %d  n_ok: %d \n", max_mds, n_ok)
	return keys // the keys sorted by value (small to large), and the number of values <= max_mds
}
