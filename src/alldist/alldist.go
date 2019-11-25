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
	"strings"
	"sort"
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
	flag.Int64Var(&seed, "seed", -1, "# rng seed (default: set from clock.)")
	flag.Float64Var(&dist_fraction, "dfrac", 1, "# fraction of all N choose 2 distances to do.")

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

	files := strings.Split(files_string, ",") // fasta files, or matrix files!

	test_file := files[0] // filename of first file, now check to see whether is fasta or matrix
	input_format := what_is_format(test_file)
	fmt.Println("# input format: ", input_format)
	if input_format == "other" {
		os.Exit(1)
	}
	//	search_time := time.Now()
	//	distance_time := search_time

	for _, file := range files {

		var sequence_set *sequenceset.Sequence_set
		if input_format == "fasta" {
			sequence_set = sequenceset.Construct_from_fasta_file(file, max_missing_data_proportion, missing_data_prob)
		} else if input_format == "matrix" {
			sequence_set = sequenceset.Construct_from_matrix_file(file, max_missing_data_proportion)
		}
		// sequence_set.Add_missing_data(missing_data_prob)
		idpair_dist := sequence_set.Fraction_of_all_distances(dist_fraction)

		

		for id1, id2_dist := range idpair_dist {
			sorted_keys := keys_sorted_by_value(id2_dist)
			fmt.Printf("%s  ", id1)
			for i2, id2 := range sorted_keys {
				d12 := id2_dist[id2]
				fmt.Printf("%s %5.3f  ", id2, d12)
				if i2 >= 19 {
					break
				}
			}
			fmt.Println("")
		}
	}
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
func keys_sorted_by_value(amap map[string]float64) ([]string) {
	keys := make([]string, 0, len(amap))
	for k, _ := range amap {
		keys = append(keys, k)
		//	fmt.Fprintf(os.Stderr, "k v: %s  %d\n", k, v)
	}
	sort.Slice(keys, func(i, j int) bool { return amap[keys[i]] < amap[keys[j]] })
//	fmt.Fprintf(os.Stderr, "max_mds %d  n_ok: %d \n", max_mds, n_ok)
	return keys // the keys sorted by value (small to large), and the number of values <= max_mds
}