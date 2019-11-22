package seqchunkset

import (
	"fmt"
	"math/rand"
	"mytypes"
	"os"
	"sequenceset"
	"sort"
	//	"testing"
)

type Sequence_chunk_set struct {
	Sequence_set              *sequenceset.Sequence_set
	Chunk_size                int
	Chunk_specs               []chunk_spec
	Chunk__seq_matchindices   map[string]map[string][]int // keys are chunk specifier strings, values are maps, whose keys are chunk seqs (e.g. '10020101') values are slices of sequence indices
	Missing_data_chunk_counts []int                       // Missing_data_chunk_counts[i] is the number of chunks with missing data for sequence with index i
	N_chunked_sequences       int                         // number of sequences entered into Chunk__seq_matchindices, Missing_data_chunk_counts so far.
}

type chunk_spec struct {
	s string   // e.g. 'A0 A34 A5 A101'
	a []string //
}

// *********************** exported functions *********************************************

func Construct_from_sequence_set(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) *Sequence_chunk_set {

	var seq_chunk_set Sequence_chunk_set
	seq_chunk_set.Chunk_specs = get_chunk_set(sequence_set, chunk_size, n_chunks)
	seq_chunk_set.Sequence_set = sequence_set
	seq_chunk_set.Chunk__seq_matchindices = make(map[string]map[string][]int)          // chunk__seq_matchindices
	seq_chunk_set.Missing_data_chunk_counts = make([]int, len(sequence_set.Sequences)) // missing_data_chunk_counts
	for i := 0; i < len(sequence_set.Sequences); i++ {
		//	id := sequence_set.SeqIndex_id[index]
		seq_chunk_set.Add_sequence()

	}
	//	}
	return &seq_chunk_set
}

// /*
func Construct_empty(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) *Sequence_chunk_set {

	var seq_chunk_set Sequence_chunk_set
	seq_chunk_set.Chunk_specs = get_chunk_set(sequence_set, chunk_size, n_chunks)
	seq_chunk_set.Sequence_set = sequence_set                                 // sequenceset.Construct_empty()
	seq_chunk_set.Chunk__seq_matchindices = make(map[string]map[string][]int) // chunk__seq_matchindices
	seq_chunk_set.Missing_data_chunk_counts = make([]int, 0)                  // missing_data_chunk_counts
	seq_chunk_set.N_chunked_sequences = 0
	return &seq_chunk_set
} /* */

func (scs *Sequence_chunk_set) Add_sequence() { // id string, sequence string) {

	new_index := scs.N_chunked_sequences // this is the number of sequences in the sequence_chunk_set so far
	seq_set := scs.Sequence_set
	sequence := seq_set.Sequences[new_index]
	//	id := seq_set.SeqIndex_id[new_index]

	scs.Missing_data_chunk_counts = append(scs.Missing_data_chunk_counts, 0)
	// new_index := len(scs.Sequence_set.Sequences) // new index is 1 greater than the max previous index; if have N sequences already, w indices 0 through N-1, new index is N

	//	scs.Sequence_set.Add_sequence(new_index, id, sequence)

	for _, chsp := range scs.Chunk_specs { // loop over chunk specifiers
		csa := chsp.a // chunk spec in array form
		css := chsp.s // chunk spec in string form
		//	css := scs.Chunk_spec_strings[ich]
		seq_matchindices, ok1 := scs.Chunk__seq_matchindices[css]
		if !ok1 {
			seq_matchindices = make(map[string][]int)
			scs.Chunk__seq_matchindices[css] = seq_matchindices
		}
		chunk_seq := get_chunk_seq(scs.Sequence_set.SnpId_index, sequence, csa, &scs.Missing_data_chunk_counts[new_index])
		matchindices, ok2 := seq_matchindices[chunk_seq] // matchindices is slice of indices of sequences which match on this chunk.
		if !ok2 {
			matchindices = make([]int, 0, 1)
			scs.Chunk__seq_matchindices[css][chunk_seq] = matchindices
		}
		scs.Chunk__seq_matchindices[css][chunk_seq] = append(matchindices, new_index) // push the sequence index onto appropriate slice
		/*  fmt.Fprintln(os.Stderr, "X")
		 for ii, mindex := range scs.Chunk__seq_matchindices[css][chunk_seq] {
			fmt.Fprintln(os.Stderr, "css:  ", css, "  chunk_seq: ", chunk_seq, "ii: ", ii, "  mindex: ", mindex)
		} /* */
	}
	scs.N_chunked_sequences++
}

// search for relative of query sequence  sequence  using the seqchunkset scs.
func (scs *Sequence_chunk_set) Get_chunk_matchindex_counts(qseq_id string, sequence string, n_top int) ([]*mytypes.IntIntIntF64, []*mytypes.IntIntIntF64, int, int) {
	seq_length := len(sequence)
	n_subj_seqs := scs.N_chunked_sequences
	n_chunks := len(scs.Chunk_specs)
	//	fmt.Fprintln(os.Stderr, "n_subj_seqs: ", n_subj_seqs, "  seq_length: ", seq_length)
	chunk_mdmd_counts := make([]int, n_subj_seqs) // chunk_mdmd_counts[i] is number of chunks in (subj) sequence i which are missing there and also in query sequence (i.e. the one whose relatives we are searching for)
	chunkwise_match_info := make([]*mytypes.IntIntIntF64, n_subj_seqs)
	for i := 0; i < n_subj_seqs; i++ {
		iiif := mytypes.IntIntIntF64{i, 0, n_chunks - scs.Missing_data_chunk_counts[i], -1}
		chunkwise_match_info[i] = &iiif
	}
	if seq_length != scs.Sequence_set.Sequence_length {
		fmt.Printf("sequence lengths: %8d %8d not equal; exiting.\n", seq_length, scs.Sequence_set.Sequence_length)
		os.Exit(1)
	}
	seq2_chunk_md_count := 0

	chunk__seq_matchindices := scs.Chunk__seq_matchindices
	for _, chsp := range scs.Chunk_specs { // loop over chunk specifiers
		csa := chsp.a
		css := chsp.s // scs.Chunk_spec_strings[ich]
		chunk_seq := ""
		for _, snp_id := range csa { // assemble the sequence chunk (string)
			idx := scs.Sequence_set.SnpId_index[snp_id]
			char := sequence[idx : idx+1]
			if char == string(mytypes.MDchar) {
				chunk_seq = "Missing_data"
				break
			}
			chunk_seq += sequence[idx : idx+1]
		}
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
						//	fmt.Fprintln(os.Stderr, "Mindex:  ", mindex)
						chunkwise_match_info[mindex].B++ // O(N^2 S_ch)
					}
				}
			}
		}
	}
	chunk_match_total_count := 0 // matches, md in neither, summed over all chunks and all subj. sequences
	chunk_mdmd_total_count := 0  //  'matches' md in both, summed over all chunks and all subj. sequences
	// A: index, B: matching chunk count, C: OK chunk count (i.e. no md), D: B/C
	chunkwise_match_info_OK := make([]*mytypes.IntIntIntF64, 0)  //
	chunkwise_match_info_BAD := make([]*mytypes.IntIntIntF64, 0) // these have
	for i, x := range chunkwise_match_info {                     // O(N^2)
		index := x.A        // index of subj. sequence
		TestEqual(i, index) // exit if not equal
		// x.B is the number of matching chunks between query and subj
		x.C -= (seq2_chunk_md_count - chunk_mdmd_counts[index]) // the number of chunks with OK data in both query and subj. seqs
		fmt.Fprintln(os.Stderr, x.C)
		if x.C <= 0 {                                           // for now, 'BAD' criterion is that there are no chunks with OK data (no md) in both query and subj.
		//	fmt.Fprintln(os.Stderr, "qseq_id: ", qseq_id, "  s index: ", index, " matchcount: ", x.B, " x.C: ", x.C)
			chunkwise_match_info_BAD = append(chunkwise_match_info_BAD, x)
			//	os.Exit(1)
		} else {
			x.D = float64(x.B) / float64(x.C) // will sort on this. (fraction of OK chunks which match)
			chunkwise_match_info_OK = append(chunkwise_match_info_OK, x)
		}
		chunk_match_total_count += x.B
		chunk_mdmd_total_count += chunk_mdmd_counts[i]
	}
	//	fmt.Println(chunk_match_total_count, chunk_mdmd_total_count)
	if n_top < len(chunkwise_match_info_OK) {
		chunkwise_match_info_OK = quickselect(chunkwise_match_info_OK, n_top) // top n_top matches, i
	}
	sort.Slice(chunkwise_match_info_OK, func(i, j int) bool { return chunkwise_match_info_OK[i].D > chunkwise_match_info_OK[j].D })
	/* for _, match_info := range chunkwise_match_info {
	   } /* */
	return chunkwise_match_info_OK, chunkwise_match_info_BAD, chunk_match_total_count, chunk_mdmd_total_count
}

// search for relatives of q_seq_set  in scs
func (scs *Sequence_chunk_set) Search(q_seq_set *sequenceset.Sequence_set, n_keep int, add bool) (map[string][]*mytypes.IntIntIntF64, map[string][]*mytypes.IntIntIntF64, int, int) {
	qid_smatchinfos := make(map[string][]*mytypes.IntIntIntF64) // keys strings (id2), values: slices
	qid_badmatches := make(map[string][]*mytypes.IntIntIntF64)
	total_chunk_match_count := 0
	total_mdmd_match_count := 0
	for qindex, qseq := range q_seq_set.Sequences {
		qid := q_seq_set.Seq_index_to_id(qindex)
		if true { // search against the previously read-in sequences
			top_smatchinfos, bad_matches, tcmc, tmdmdc := scs.Get_chunk_matchindex_counts(qid, qseq, n_keep)
			total_chunk_match_count += tcmc
			total_mdmd_match_count += tmdmdc
			if qindex%1000 == 0 {
				fmt.Fprintf(os.Stderr, "Search %d done.\n", qindex)
			}
			qid_smatchinfos[qid] = top_smatchinfos
			if len(bad_matches) > 0 {
				qid_badmatches[qid] = bad_matches
			}
		}
		if add {
			scs.Add_sequence() // add latest sequence
		}
	}
	return qid_smatchinfos, qid_badmatches, total_chunk_match_count, total_mdmd_match_count
}

// ***************************************************************************************************

// ******************* functions which are not exported **********************************************

// get a set of n_chunks each of size chunk_size
func get_chunk_set(sequence_set *sequenceset.Sequence_set, chunk_size int, n_chunks int) []chunk_spec {

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	// and also the corresponding chunk-specifying slices (e.g. []int{1,5,109,4} )

	chunk_specs := make([]chunk_spec, 0)
	if n_chunks < 0 {
		n_chunks = sequence_set.N_ok_snps / chunk_size //n_chunks_from_sequence // just do one pass through the set of sequences
	}
	fmt.Fprintf(os.Stderr, "N_ok_snps: %d %d\n", sequence_set.N_ok_snps, n_chunks)
	for k := 0; true; k++ {
		is := make([]int, sequence_set.N_ok_snps) // Sequence_length) // slice of integers
		for j, _ := range is {
			is[j] = int(j)
		}
		if k >= 0 {
			rand.Shuffle(len(is), func(i, j int) { is[i], is[j] = is[j], is[i] })
		}
		// is is now a slice of sequence_length integers, (0 through sequence_length-1) in randomized order.
		for used_so_far := 0; used_so_far+chunk_size <= sequence_set.N_ok_snps; used_so_far += chunk_size {
			//	fmt.Fprintln(os.Stderr, used_so_far, is[used_so_far])

			snp_id := sequence_set.Sorted_snp_ids[is[used_so_far]]
			//	fmt.Fprintln(os.Stderr, snp_id, used_so_far, is[used_so_far])
			//	fmt.Fprintf(os.Stderr, "xxx %s ", snp_id)
			chunk_spec_str := snp_id // fmt.Sprintf("_%d", next_gt)

			var chunk_spec_array []string                       // int
			chunk_spec_array = append(chunk_spec_array, snp_id) // next_gt)  // snp_id)
			for i := 1; i < chunk_size; i++ {
				//		next_gt = is[used_so_far + i]
				snp_id = sequence_set.Sorted_snp_ids[is[used_so_far+i]]
				//	fmt.Fprintf(os.Stderr, "%s ", snp_id)

				chunk_spec_str += " " + snp_id                      // fmt.Sprintf("_%d", next_gt)
				chunk_spec_array = append(chunk_spec_array, snp_id) // next_gt)
			}
			//	fmt.Fprintf(os.Stderr, "\n")
			chunk_specs = append(chunk_specs, chunk_spec{chunk_spec_str, chunk_spec_array})
			if len(chunk_specs) >= n_chunks {
				return chunk_specs
			}
		}
	}
	os.Exit(1)
	return chunk_specs
}

func get_chunk_seq(snpid_index map[string]int, seq string, csa []string, md_count *int) string {
	chunk_seq := ""
	for _, snp_id := range csa {
		snp_index := snpid_index[snp_id]
		char := seq[snp_index : snp_index+1]
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
			px := len(list) / 2                         // 'middle' index
			pv := list[px].D                            // 'pivot' value
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

func TestEqual(a int, b int) {

	if a != b {
		os.Exit(1)
	}
}


// incrementally search and construct the seqchset, i.e. search each seq. against existing seqchset, then
// add it to the seqchset.
/* func (scs *Sequence_chunk_set) Search_and_construct(n_keep int) (map[string][]*mytypes.IntIntIntF64, map[string][]*mytypes.IntIntIntF64, int, int) {
	qid_smatchinfos := make(map[string][]*mytypes.IntIntIntF64) // keys strings (id2), values: slices
	qid_badmatches := make(map[string][]*mytypes.IntIntIntF64)
	total_chunk_match_count := 0
	total_mdmd_match_count := 0
	for qindex, qseq := range scs.Sequence_set.Sequences {
		qid := scs.Sequence_set.Seq_index_to_id(qindex)
		if true { // search against the previously read-in sequences
			top_smatchinfos, bad_matches, tcmc, tmdmdc := scs.Get_chunk_matchindex_counts(qid, qseq, n_keep)
			total_chunk_match_count += tcmc
			total_mdmd_match_count += tmdmdc
			if qindex%1000 == 0 {
				fmt.Fprintf(os.Stderr, "Search %d done.\n", qindex)
			}
			qid_smatchinfos[qid] = top_smatchinfos
			if len(bad_matches) > 0 {
				qid_badmatches[qid] = bad_matches
			}
		}
		scs.Add_sequence() // add latest sequence
	}
	return qid_smatchinfos, qid_badmatches, total_chunk_match_count, total_mdmd_match_count
} /* */
