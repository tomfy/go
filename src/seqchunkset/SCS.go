package seqchunkset

import (
	//	"bufio"
	"fmt"
	"os"
	//	"reflect"
	//	"regexp"
	"sequenceset"
	"sort"
	//	"strings"
)

type sequence_chunk_set struct {
	//	Sequences               []string
	//	Sequence_length         int
	//	Index_id                map[int]string
	//	Id_index                map[string]int
	Sequence_set            *sequenceset.Sequence_set
	Chunk_size              int
	Chunk_spec_strings      []string
	Chunk_spec_arrays       [][]int
	Chunk__seq_matchindices map[string]map[string][]int // keys are chunk specifier strings, values are pts to maps, whose keys are chunk seqs (e.g. '10020101') values are
}

func Construct_from_fasta_file(filename string, chunk_size int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set
	chunk__seq_matchindices := make(map[string]map[string][]int)
	//	var sequences []string
	var sequence_set *sequenceset.Sequence_set

	sequence_set = sequenceset.Construct_from_fasta_file(filename)

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	var chunk_spec_strings []string
	var chunk_spec_arrays [][]int
	for used_so_far := 0; used_so_far+chunk_size <= sequence_set.Sequence_length; used_so_far += chunk_size {
		chunk_spec_str := fmt.Sprintf("%d", used_so_far)
		var chunk_spec_array []int
		chunk_spec_array = append(chunk_spec_array, used_so_far)
		for i := 1; i < chunk_size; i++ {
			chunk_spec_str += fmt.Sprintf("_%d", used_so_far+i)
			chunk_spec_array = append(chunk_spec_array, used_so_far+i)
		}
		chunk_spec_strings = append(chunk_spec_strings, chunk_spec_str)
		chunk_spec_arrays = append(chunk_spec_arrays, chunk_spec_array)
	}
	seq_chunk_set.Chunk_spec_strings = chunk_spec_strings
	seq_chunk_set.Chunk_spec_arrays = chunk_spec_arrays

	// for each sequence, consider as chunks and store sequence index in
	for ich, csa := range chunk_spec_arrays { // loop over chunk specifiers
		css := chunk_spec_strings[ich]
		seq_matchindices, ok1 := chunk__seq_matchindices[css]
		//	fmt.Println(reflect.TypeOf(seq_matchindices))
		if !ok1 {
			seq_matchindices = make(map[string][]int)
			chunk__seq_matchindices[css] = seq_matchindices
		}
		for iseq, seq := range sequence_set.Sequences { // loop over the sequences
			//	seq := sequences[iseq]
			//	seq_as_array := strings.Split(seq, "")
			//	fmt.Println("XXX: " + seq)
			chunk_seq := ""
			for j := 0; j < chunk_size; j++ { // assemble the sequence chunk
				//	snp := seq[csa[j] : csa[j]+1]

				chunk_seq += seq[csa[j] : csa[j]+1]
				//	fmt.Printf("AAAAA: %d   %s\n", j, chunk_seq)
			}
			//	Chunk__seq_matchindices *map[string]*map[string]*[]int
			//	chunk__seq_matchindices := make(map[string]*map[string]*[]int)
			matchindices, ok2 := seq_matchindices[chunk_seq]
			if !ok2 {
				matchindices = make([]int, 0, 1)
				chunk__seq_matchindices[css][chunk_seq] = matchindices
			}
			chunk__seq_matchindices[css][chunk_seq] = append(matchindices, iseq) // push the sequence index onto appropriate slice
		}
	}

	seq_chunk_set.Sequence_set = sequence_set
	seq_chunk_set.Chunk__seq_matchindices = chunk__seq_matchindices
	return &seq_chunk_set
}

func (scs sequence_chunk_set) Get_chunk_matchindex_counts(sequence string, matchindex_counts []int) ([]int, []int) {
	seq_length := len(sequence)
	if seq_length != scs.Sequence_set.Sequence_length {
		fmt.Printf("sequence lengths: %8d %8d not equal; exiting.\n", seq_length, scs.Sequence_set.Sequence_length)
		os.Exit(1)
	}
	//	fmt.Println("# "+sequence)
	chunk__seq_matchindices := scs.Chunk__seq_matchindices
	for ich, csa := range scs.Chunk_spec_arrays { // loop over chunk specifiers
		css := scs.Chunk_spec_strings[ich]
		//fmt.Printf("%d  %s \n", ich, css)
		//fmt.Println(csa)
		chunk_seq := ""
		for _, idx := range csa { // assemble the sequence chunk (string)
			//	fmt.Printf("GGG:  %s %d\n", css, idx)
			chunk_seq += sequence[idx : idx+1]
			//		fmt.Printf("HHH:  %s %d\n", chunk_seq, idx)
		}
		//	fmt.Println("JJJJJ")
	//	fmt.Printf("#   %s  %s\n", css, chunk_seq)
		seq_matchindices, ok1 := chunk__seq_matchindices[css]
		if ok1 {
			matchindices, ok2 := seq_matchindices[chunk_seq]
			if ok2 {
				for _, mindex := range matchindices {
				//	fmt.Printf("      %d  %d\n", j, mindex)
					matchindex_counts[mindex]++
				}
			}
		}
	}

	fmt.Printf("XXX:   ")
	for midx, counts := range matchindex_counts{
		fmt.Printf("%d %d   ", midx, counts)
	}
	fmt.Printf("\n\n")

	n_seqs := len(matchindex_counts)
	//	sort.Ints(matchindex_counts)                 // small to large
	indices := make([]int, n_seqs)
	for i, _ := range indices {
		indices[i] = i
	}
	sort.Slice(indices, func(i, j int) bool { return matchindex_counts[i] > matchindex_counts[j] })
	counts := make([]int, 15)
	for i, idx := range indices[0:15] {
		counts[i] = matchindex_counts[idx]
		fmt.Printf("%2d %2d %2d  ", i, idx, counts[i])
	}
	fmt.Printf("\n")
	return indices[0:15], counts
}
