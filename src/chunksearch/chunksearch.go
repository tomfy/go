package main

import (
//	"bufio"
	"flag"
	"fmt"
//	"os"
//	"regexp"
//	"strings"
	"seqchunkset"
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

	file1Ptr := flag.String("f1", "", "name of first fasta file.")
	file2Ptr := flag.String("f2", "", "name of 2nd fasta file.")
	chunk_size_ptr := flag.Int("size", 6, "number of snps in each chunk")

	flag.Parse()

	fmt.Printf("%s\n", *file1Ptr)
	fmt.Printf("%s\n", *file2Ptr)

	seqchset := seqchunkset.Construct_from_fasta_file(*file1Ptr, *chunk_size_ptr)
	//	fmt.Printf("%3d\n", *f1retval)
	fmt.Println(seqchset.Sequences)
	fmt.Println(seqchset.Index_id)
	fmt.Println(seqchset.Id_index)


	fmt.Println(seqchset.Chunk_spec_strings)

		
	/*
	   	fh2, err := os.Open(*file2Ptr)
	   	if err != nil {
	   	os.Exit(1)
	   }
	*/

}
