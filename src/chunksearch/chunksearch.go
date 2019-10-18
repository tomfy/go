package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"regexp"
	"strings"
	"fasta"
)

func store_fasta_file(filename string) int {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	
	var sequences1 []string
	index_id := make(map[int]string)
	id_index := make(map[string]int)
	/*	var chunk_specifier_strings []string  /* elements like '23_51_3_1_135' */
	min_seq_len := 1000000000
	max_seq_len := -1

	index := 0
	scanner1 := bufio.NewScanner(fh)
	for scanner1.Scan() {
		line := scanner1.Text()
		r, _ := regexp.Compile("^>([0-9]+).*")
		if r.MatchString(line) { /* >id line */
			match_strings := r.FindStringSubmatch(line)
			index_id[index] = match_strings[1]
			id_index[match_strings[1]] = index
			fmt.Println(match_strings[1])
		} else { /* sequence line */
			line = strings.TrimSpace(line)
			if len(line) < min_seq_len {
				min_seq_len = len(line)
			}
			if len(line) > max_seq_len {
				max_seq_len = len(line)
			}
			sequences1 = append(sequences1, line)
			fmt.Printf("   [%s]\n", line)
			index++
		}
	}
	if min_seq_len != max_seq_len {
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}
	return 1
}

func main() {

	/* command line options: */

	/* input files: */
	file1Ptr := flag.String("f1", "", "name of first fasta file.")
	file2Ptr := flag.String("f2", "", "name of 2nd fasta file.")

	flag.Parse()

	fmt.Printf("%s\n", *file1Ptr)
	fmt.Printf("%s\n", *file2Ptr)

	f1retval := store_fasta_file(*file1Ptr)
	fmt.Printf("%d8\n", f1retval)

/*	var sequences1 []string
	index_id := make(map[int]string)
	id_index := make(map[string]int)
	//	var chunk_specifier_strings []string  // elements like '23_51_3_1_135' 
	min_seq_len := 1000000000
	max_seq_len := -1

	index := 0
	scanner1 := bufio.NewScanner(fh1)
	for scanner1.Scan() {
		line := scanner1.Text()
		r, _ := regexp.Compile("^>([0-9]+).*")
		if r.MatchString(line) { // >id line 
			match_strings := r.FindStringSubmatch(line)
			index_id[index] = match_strings[1]
			id_index[match_strings[1]] = index
			fmt.Println(match_strings[1])
		} else { // sequence line 
			line = strings.TrimSpace(line)
			if len(line) < min_seq_len {
				min_seq_len = len(line)
			}
			if len(line) > max_seq_len {
				max_seq_len = len(line)
			}
			sequences1 = append(sequences1, line)
			fmt.Printf("   [%s]\n", line)
			index++
		}
	}
	if min_seq_len != max_seq_len {
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}

	*/
		
	/*
	   	fh2, err := os.Open(*file2Ptr)
	   	if err != nil {
	   	os.Exit(1)
	   }
	*/

}
