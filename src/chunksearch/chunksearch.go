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

/*
	if min_seq_len != max_seq_len {
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}
	return 1
} */

func main() {

	/* command line options: */

	/* input files: */

	file1 := flag.String("f1", "", "name of first fasta file.")
	file2 := flag.String("f2", "", "name of 2nd fasta file.")
	chunk_size_ptr := flag.Int("size", 6, "number of snps in each chunk")

	flag.Parse()

	fmt.Printf("file1: %s\n", *file1)
	fmt.Printf("file2: %s\n", *file2)

	seqchset := seqchunkset.Construct_from_fasta_file(*file1, *chunk_size_ptr)
	//	fmt.Printf("%3d\n", *f1retval)
	fmt.Println(seqchset.Sequence_set.Sequences)
	fmt.Println(seqchset.Sequence_set.Index_id)
	fmt.Println(seqchset.Sequence_set.Id_index)
	fmt.Println(seqchset.Chunk_spec_strings)

	for i, css := range seqchset.Chunk_spec_strings {
fmt.Printf("TTT: %d\n", i)
		seq_matchindices := seqchset.Chunk__seq_matchindices[css]
		for chunkseq, _ := range seq_matchindices {
			fmt.Printf("SSS: %d   %s   %s\n", i, css, chunkseq )
		}
	}

	sequence_set2 := sequenceset.Construct_from_fasta_file(*file2)
	fmt.Println("WWWW")
	for _, seq2 := range sequence_set2.Sequences {
		matchindex_counts := make([]int, 1, 1)
		fmt.Printf("%d %d \n", len(matchindex_counts), cap(matchindex_counts))
		highest_counts := seqchset.Get_chunk_matchindex_counts(seq2, matchindex_counts)
		fmt.Println(highest_counts)
	}
	/*
	   	fh2, err := os.Open(*file2Ptr)
	   	if err != nil {
	   	os.Exit(1)
	   }
	*/

}
