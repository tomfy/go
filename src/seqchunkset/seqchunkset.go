package seqchunkset

import (
	"fmt"
	"mytypes"
	"os"
	"sequenceset"
	"sort"
)

type sequence_chunk_set struct {
	Sequence_set            *sequenceset.Sequence_set
	Chunk_size              int
	Chunk_spec_strings      []string
	Chunk_spec_arrays       [][]int
	Chunk__seq_matchindices map[string]map[string][]int // keys are chunk specifier strings, values are pts to maps, whose keys are chunk seqs (e.g. '10020101') values are
}

func Construct_from_fasta_file(sequence_set *sequenceset.Sequence_set, chunk_size int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set
	chunk__seq_matchindices := make(map[string]map[string][]int)

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	// and also the corresponding chunk-specifying slices (e.g. []int{1,5,109,4} )
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
		if !ok1 {
			seq_matchindices = make(map[string][]int)
			chunk__seq_matchindices[css] = seq_matchindices
		}
		for iseq, seq := range sequence_set.Sequences { // loop over the sequences

			chunk_seq := ""
			for j := 0; j < chunk_size; j++ { // assemble the sequence chunk
				chunk_seq += seq[csa[j] : csa[j]+1]
			}
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

func (scs sequence_chunk_set) Get_chunk_matchindex_counts(sequence string, mindex_count_pairs []mytypes.Pair_int_int, n_top int) []mytypes.Pair_int_int {

//	fmt.Println(mindex_count_pairs)
	seq_length := len(sequence)
	if seq_length != scs.Sequence_set.Sequence_length {
		fmt.Printf("sequence lengths: %8d %8d not equal; exiting.\n", seq_length, scs.Sequence_set.Sequence_length)
		os.Exit(1)
	}
	chunk__seq_matchindices := scs.Chunk__seq_matchindices
	for ich, csa := range scs.Chunk_spec_arrays { // loop over chunk specifiers
		css := scs.Chunk_spec_strings[ich]
		chunk_seq := ""
		for _, idx := range csa { // assemble the sequence chunk (string)
			chunk_seq += sequence[idx : idx+1]
		}
		seq_matchindices, ok1 := chunk__seq_matchindices[css]
		if ok1 {
			matchindices, ok2 := seq_matchindices[chunk_seq]
			if ok2 {
				for _, mindex := range matchindices {
					mindex_count_pairs[mindex].B++
				}
			}
		}
	}
	sort.Slice(mindex_count_pairs, func(i, j int) bool { return mindex_count_pairs[i].B > mindex_count_pairs[j].B })
	return mindex_count_pairs[0:n_top-1]
}
