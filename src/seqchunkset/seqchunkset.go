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
	Sequence_set              *sequenceset.Sequence_set
	Chunk_size                int
	Chunk_spec_strings        []string
	Chunk_spec_arrays         [][]int
	Chunk__seq_matchindices   map[string]map[string][]int // keys are chunk specifier strings, values are maps, whose keys are chunk seqs (e.g. '10020101') values are slices of sequence indices
	Missing_data_chunk_counts []int                       // Missing_data_chunk_counts[i] is the number of chunks with missing data for sequence with index i
}

// *********************** exported functions *********************************************

func Construct_from_sequence_set(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set

	chunk_spec_strings, chunk_spec_arrays := get_chunk_set(sequence_set.Sequence_length, chunk_size, n_chunks)
	seq_chunk_set.Chunk_spec_strings = chunk_spec_strings
	seq_chunk_set.Chunk_spec_arrays = chunk_spec_arrays

	if false {
		chunk__seq_matchindices := make(map[string]map[string][]int)
		missing_data_chunk_counts := make([]int, len(sequence_set.Sequences))
		// for each sequence, consider as chunks and store sequence index in
		for ich, csa := range chunk_spec_arrays { // loop over chunk specifiers
			css := chunk_spec_strings[ich]
			seq_matchindices, ok1 := chunk__seq_matchindices[css]
			if !ok1 {
				seq_matchindices = make(map[string][]int)
				chunk__seq_matchindices[css] = seq_matchindices
			}
			for iseq, seq := range sequence_set.Sequences { // loop over the sequences
				chunk_seq := get_chunk_seq(seq, csa, &missing_data_chunk_counts[iseq])
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
		seq_chunk_set.Missing_data_chunk_counts = missing_data_chunk_counts
	} else {
		// index_offset := len(sequence_set.Sequences)
		seq_chunk_set.Sequence_set = sequence_set
		seq_chunk_set.Chunk__seq_matchindices = make(map[string]map[string][]int)          // chunk__seq_matchindices
		seq_chunk_set.Missing_data_chunk_counts = make([]int, len(sequence_set.Sequences)) // missing_data_chunk_counts
		for index, seq := range sequence_set.Sequences {
			id := sequence_set.Index_id[index]
			seq_chunk_set.Add_sequence(id, seq)
		}
	}
	return &seq_chunk_set
}

func Construct_empty(sequence_length int, chunk_size int, n_chunks int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set

	chunk_spec_strings, chunk_spec_arrays := get_chunk_set(sequence_length, chunk_size, n_chunks)
	seq_chunk_set.Chunk_spec_strings = chunk_spec_strings
	seq_chunk_set.Chunk_spec_arrays = chunk_spec_arrays

	seq_chunk_set.Sequence_set = sequenceset.Construct_empty()
	seq_chunk_set.Chunk__seq_matchindices = make(map[string]map[string][]int)                        // chunk__seq_matchindices
	seq_chunk_set.Missing_data_chunk_counts = make([]int, len(seq_chunk_set.Sequence_set.Sequences)) // missing_data_chunk_counts
	return &seq_chunk_set
}

//*
func (scs *sequence_chunk_set) Add_sequence(id string, sequence string) {
	//	chunk_spec_strings := scs.Chunk_spec_strings
	//	chunk_spec_arrays := scs.Chunk_spec_arrays
	//	chunk__seq_matchindices := scs.Chunk__seq_matchindices
	//	missing_data_chunk_counts := scs.Missing_data_chunk_counts
	scs.Missing_data_chunk_counts = append(scs.Missing_data_chunk_counts, 0)
	new_index := len(scs.Sequence_set.Sequences) // new index is 1 greater than the max previous index; if have N sequences already, w indices 0 through N-1, new index is N
	//	fmt.Printf("xxx %d \n", len(scs.Sequence_set.Sequences))
	scs.Sequence_set.Add_sequence(new_index, id, sequence)
	//	fmt.Printf("%d %s %d %s\n", new_index, id, scs.Sequence_set.Id_index[id], scs.Sequence_set.Index_id[new_index])

	for ich, csa := range scs.Chunk_spec_arrays { // loop over chunk specifiers
		css := scs.Chunk_spec_strings[ich]
		seq_matchindices, ok1 := scs.Chunk__seq_matchindices[css]
		if !ok1 {
			seq_matchindices = make(map[string][]int)
			scs.Chunk__seq_matchindices[css] = seq_matchindices
		}
		//	for iseq, seq := range sequence_set.Sequences { // loop over the sequences

		chunk_seq := get_chunk_seq(sequence, csa, &scs.Missing_data_chunk_counts[new_index])
		matchindices, ok2 := seq_matchindices[chunk_seq] // matchindices is slice of indices of sequences which match on this chunk.
		if !ok2 {
			matchindices = make([]int, 0, 1)
			scs.Chunk__seq_matchindices[css][chunk_seq] = matchindices
		}
		scs.Chunk__seq_matchindices[css][chunk_seq] = append(matchindices, new_index) // push the sequence index onto appropriate slice
		//	}
	}

} /* */

func (scs *sequence_chunk_set) Get_chunk_matchindex_counts(sequence string /* chunkwise_match_info []mytypes.IntIntIntF64,*/, n_top int) []*mytypes.IntIntIntF64 {

	//	fmt.Println(chunkwise_match_info)
	seq_length := len(sequence)

	n_seqs := len(scs.Sequence_set.Sequences)
	n_chunks := len(scs.Chunk_spec_strings)

	chunkwise_match_info := make([]*mytypes.IntIntIntF64, n_seqs)
	for i := 0; i < n_seqs; i++ {
		iiif := mytypes.IntIntIntF64{i, 0, n_chunks - scs.Missing_data_chunk_counts[i], -1}
		chunkwise_match_info[i] = &iiif
		//	chunkwise_match_info[i].A = i
		//	chunkwise_match_info[i].C = n_chunks - scs.Missing_data_chunk_counts[i] // number of OK chunks in seq i (OK == no missing data)
	}
	if seq_length != scs.Sequence_set.Sequence_length {
		fmt.Printf("sequence lengths: %8d %8d not equal; exiting.\n", seq_length, scs.Sequence_set.Sequence_length)
		os.Exit(1)
	}
	seq2_chunk_md_count := 0
	chunk_mdmd_counts := make([]int, n_seqs) // chunk_mdmd_counts[i] is number of chunks in (db) sequence i which are missing there and also in sequence (i.e. the one whose relatives we are searching for)
	chunk__seq_matchindices := scs.Chunk__seq_matchindices
	for ich, csa := range scs.Chunk_spec_arrays { // loop over chunk specifiers
		css := scs.Chunk_spec_strings[ich]
		chunk_seq := ""
		for _, idx := range csa { // assemble the sequence chunk (string)
			char := sequence[idx : idx+1]
			if char == string(mytypes.MDchar) {
				chunk_seq = "Missing_data"
				break
			}
			chunk_seq += sequence[idx : idx+1]
		}
		//	fmt.Println(chunk_seq)
		seq_matchindices, ok1 := chunk__seq_matchindices[css]
		if ok1 {
			matchindices, ok2 := seq_matchindices[chunk_seq]
			if ok2 {
				if chunk_seq == "Missing_data" {
					seq2_chunk_md_count++ // count the chunks with missing data in sequence
					for _, mindex := range matchindices {
						chunk_mdmd_counts[mindex]++
					}
				} else {
					for _, mindex := range matchindices {
						chunkwise_match_info[mindex].B++
					}
				}
			}
		}
	}
	for _, x := range chunkwise_match_info {
		index := x.A
		x.C -= (seq2_chunk_md_count - chunk_mdmd_counts[index]) // this is now the number of chunks with OK data in both query and subj. sequences
		x.D = float64(x.B) / float64(x.C)                       // will sort on this. ( fraction of OK chunks which match )
		//	fmt.Printf("# %d %d %d %g   %d %d \n", x.A, x.B, chunkwise_match_info[i].C, chunkwise_match_info[i].D, seq2_chunk_md_count, chunk_mdmd_counts[index])
	}

	if n_top < n_seqs {
		chunkwise_match_info = quickselect(chunkwise_match_info, n_top) // top n_top matches, i
	}
	sort.Slice(chunkwise_match_info, func(i, j int) bool { return chunkwise_match_info[i].D > chunkwise_match_info[j].D })
	/* for _, xxx := range chunkwise_match_info {
		fmt.Println(*xxx)
	} */
	return chunkwise_match_info
}

// ***************************************************************************************************

// ******************* functions which are not exported **********************************************

// get a set of n_chunks each of size chunk_size
func get_chunk_set(sequence_length int, chunk_size int, n_chunks int) ([]string, [][]int) {

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	// and also the corresponding chunk-specifying slices (e.g. []int{1,5,109,4} )
	var chunk_spec_strings []string
	var chunk_spec_arrays [][]int
	if n_chunks < 0 {
		n_chunks = sequence_length / chunk_size //n_chunks_from_sequence // just do one pass through the set of sequences
	}
	//	fmt.Fprintf(os.Stderr, "# n_chunks: %d\n", n_chunks)

	for k := 0; true; k++ {
		is := make([]int, sequence_length) // slice of integers
		for j, _ := range is {
			is[j] = j
		}
		if k >= 0 {
			rand.Shuffle(len(is), func(i, j int) { is[i], is[j] = is[j], is[i] })
		}
		for used_so_far := 0; used_so_far+chunk_size <= sequence_length; used_so_far += chunk_size {
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

func get_chunk_seq(seq string, csa []int, md_count *int) string {
	chunk_seq := ""
	for j, _ := range csa {
		//		chunk_chars[j] = chars[csa[j]]
		char := seq[csa[j] : csa[j]+1]
		if char == string(mytypes.MDchar) {
			chunk_seq = "Missing_data"
			*md_count++
			break
		} else {
			chunk_seq += char
		}
	}
	return chunk_seq
}

// find the best k, i.e. the k with list[i].D the largest
func quickselect(list []*mytypes.IntIntIntF64, k int) []*mytypes.IntIntIntF64 {
	if k > 0 {
		k--
		keepers := make([]*mytypes.IntIntIntF64, 0, k)
		for {
			// partition
			px := len(list) / 2 // 'middle' index
			// pv := list[px].B                           // 'pivot' value
			pv := list[px].D
			last := len(list) - 1                       // index of last element
			list[px], list[last] = list[last], list[px] // swap 'middle' and last elements
			i := 0
			for j := 0; j < last; j++ { // for each index up to but not including the last ...
				if list[j].D > pv { // if elem j is less than the pivot value ...
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
		return make([]*mytypes.IntIntIntF64, 0, 0)
	}
}
