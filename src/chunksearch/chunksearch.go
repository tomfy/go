package main

import (
	//	"bufio"
	"flag"
	"fmt"
	//	"os"
	//	"regexp"
	//	"strings"
	"seqchunkset"
	"sequenceset"
)

func main() {

	/* command line options: */

	/* input files: */

	file1 := flag.String("f1", "", "name of first fasta file.")
	file2 := flag.String("f2", "", "name of 2nd fasta file.")
	chunk_size_ptr := flag.Int("size", 6, "number of snps in each chunk")

	flag.Parse()

	fmt.Printf("# file1: %s\n", *file1)
	fmt.Printf("# file2: %s\n", *file2)

	seqchset := seqchunkset.Construct_from_fasta_file(*file1, *chunk_size_ptr)
	n_seqs1 := len(seqchset.Sequence_set.Sequences)


	sequence_set2 := sequenceset.Construct_from_fasta_file(*file2)
	for index, seq2 := range sequence_set2.Sequences {
		id2 := sequence_set2.Index_to_id(index)
		fmt.Printf("%8s  ", id2)
		matchindex_counts := make([]int, n_seqs1)
		//	fmt.Printf("%d %d \n", len(matchindex_counts), cap(matchindex_counts))
		best_match_indices, counts := seqchset.Get_chunk_matchindex_counts(seq2, matchindex_counts)
		//	fmt.Printf("# ")
	//	for i, index1 := range best_match_indices {
		//	id1 := sequence_set1.Index_to_id(index1)
	//		fmt.Printf("%2d %2d  ", index1, counts[i])
	//	}
fmt.Println(best_match_indices)
fmt.Println(counts)
//		fmt.Printf("\n")
		//fmt.Println(highest_counts)
	}

}
