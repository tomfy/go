package seqchunkset

import (
	"fmt"
	"math/rand"
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

func Construct_from_sequence_set(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set
	chunk__seq_matchindices := make(map[string]map[string][]int)

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	// and also the corresponding chunk-specifying slices (e.g. []int{1,5,109,4} )
	//	var chunk_spec_strings []string
	//	var chunk_spec_arrays [][]int
	//	var n_full_chunk_sets int
	/*	seq_length := sequence_set.Sequence_length
		n_chunks_from_sequence := seq_length / chunk_size // int division here
		if n_chunks < 0 {
			n_chunks = n_chunks_from_sequence // just do one pass through the set of sequences
			n_full_chunk_sets = 1
		} else {
			n_full_chunk_sets = n_chunks / n_chunks_from_sequence
		}
		fmt.Fprintf(os.Stderr, "# full chunk sets: %d\n", n_full_chunk_sets)
		for k := 0; k < n_full_chunk_sets; k++ {
			is := make([]int, seq_length, seq_length) // slice of integers
			for j, _ := range is {
				is[j] = j
			}
			if k > 0 {
				rand.Shuffle(len(is), func(i, j int) { is[i], is[j] = is[j], is[i] })
			}
			for used_so_far := 0; used_so_far+chunk_size <= sequence_set.Sequence_length; used_so_far += chunk_size {
				next_gt := is[used_so_far]
				chunk_spec_str := fmt.Sprintf("%d", next_gt)
				var chunk_spec_array []int
				chunk_spec_array = append(chunk_spec_array, next_gt)
				for i := 1; i < chunk_size; i++ {
					next_gt = used_so_far + i
					chunk_spec_str += fmt.Sprintf("_%d", next_gt)
					chunk_spec_array = append(chunk_spec_array, next_gt)
				}
				chunk_spec_strings = append(chunk_spec_strings, chunk_spec_str)
				chunk_spec_arrays = append(chunk_spec_arrays, chunk_spec_array)

			}
		}
		seq_chunk_set.Chunk_spec_strings = chunk_spec_strings
		seq_chunk_set.Chunk_spec_arrays = chunk_spec_arrays */

	chunk_spec_strings, chunk_spec_arrays := Get_chunk_set(sequence_set, chunk_size, n_chunks)
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

func Get_chunk_set(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) ([]string, [][]int) {
	//	var seq_chunk_set sequence_chunk_set
	//	chunk__seq_matchindices := make(map[string]map[string][]int)

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	// and also the corresponding chunk-specifying slices (e.g. []int{1,5,109,4} )
	var chunk_spec_strings []string
	var chunk_spec_arrays [][]int
//	var n_full_chunk_sets int
	seq_length := sequence_set.Sequence_length
	if n_chunks < 0 {
		n_chunks = seq_length / chunk_size //n_chunks_from_sequence // just do one pass through the set of sequences
	}
	fmt.Fprintf(os.Stderr, "# n_chunks: %d\n", n_chunks)

	for k := 0; true; k++ {
		is := make([]int, seq_length, seq_length) // slice of integers
		for j, _ := range is {
			is[j] = j
		}
		if k > 0 {
			rand.Shuffle(len(is), func(i, j int) { is[i], is[j] = is[j], is[i] })
		}
		for used_so_far := 0; used_so_far+chunk_size <= sequence_set.Sequence_length; used_so_far += chunk_size {
			next_gt := is[used_so_far]
			chunk_spec_str := fmt.Sprintf("%d", next_gt)
			var chunk_spec_array []int
			chunk_spec_array = append(chunk_spec_array, next_gt)
			for i := 1; i < chunk_size; i++ {
				next_gt = used_so_far + i
				chunk_spec_str += fmt.Sprintf("_%d", next_gt)
				chunk_spec_array = append(chunk_spec_array, next_gt)
			}
			chunk_spec_strings = append(chunk_spec_strings, chunk_spec_str)
			chunk_spec_arrays = append(chunk_spec_arrays, chunk_spec_array)
			if len(chunk_spec_strings) >= n_chunks {
				return chunk_spec_strings, chunk_spec_arrays
			}
		}
	}
	os.Exit(1)
	return chunk_spec_strings, chunk_spec_arrays
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
	if true {
		mindex_count_pairs = quickselect(mindex_count_pairs, n_top)
		sort.Slice(mindex_count_pairs, func(i, j int) bool { return mindex_count_pairs[i].B > mindex_count_pairs[j].B })
		return mindex_count_pairs
	} else {
		sort.Slice(mindex_count_pairs, func(i, j int) bool { return mindex_count_pairs[i].B > mindex_count_pairs[j].B })
		return mindex_count_pairs[0 : n_top-1]
	}
}

func quickselect(list []mytypes.Pair_int_int, k int) []mytypes.Pair_int_int {
	if k > 0 {
		k--
		keepers := make([]mytypes.Pair_int_int, 0, k)
		for {
			// partition
			px := len(list) / 2                         // 'middle' index
			pv := list[px].B                            // 'pivot' value
			last := len(list) - 1                       // index of last element
			list[px], list[last] = list[last], list[px] // swap 'middle' and last elements
			i := 0
			for j := 0; j < last; j++ { // for each index up to but not including the last ...
				if list[j].B > pv { // if elem j is less than the pivot value ...
					list[i], list[j] = list[j], list[i] // swap elem j to ith position (i <= j)
					i++
				}
			} // now all elements < pv (i of them) are on left, >= pv on right
			// select
			if k == i {
				list[i], list[last] = list[last], list[i]
				keepers = append(keepers, list[:i+1]...)
				return keepers
			}
			if k < i {
				list = list[:i] // kth smallest will be one of the first i
			} else { // k > i;  kth smallest will be one of the k-(i+1) elements >= pv
				list[i], list[last] = list[last], list[i]
				keepers = append(keepers, list[:i+1]...) // 0 through ith elements are keepers
				list = list[i+1:]
				k -= i + 1
			}
		}
	} else {
		return make([]mytypes.Pair_int_int, 0, 0)
	}
}
