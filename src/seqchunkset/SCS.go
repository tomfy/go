package seqchunkset

import (
	"bufio"
	"fmt"
	"os"
	"reflect"
	"regexp"
	"sort"
	"strings"
	"sequenceset"
)

type sequence_chunk_set struct {
	Sequences       []string
	Sequence_length int
	Index_id        map[int]string
	Id_index        map[string]int
	Sequence_set sequenceset.sequence_set
	Chunk_size              int
	Chunk_spec_strings      []string
	Chunk_spec_arrays       [][]int
	Chunk__seq_matchindices map[string]map[string][]int // keys are chunk specifier strings, values are pts to maps, whose keys are chunk seqs (e.g. '10020101') values are
}

func Construct_from_fasta_file(filename string, chunk_size int) *sequence_chunk_set {

	var seq_chunk_set sequence_chunk_set
	chunk__seq_matchindices := make( map[string]map[string][]int )
	var sequences []string
	var sequence_set sequenceset
	if true {
		fh, err := os.Open(filename)
		if err != nil {
			os.Exit(1)
		}

		id_index := make(map[string]int)
		index_id := make(map[int]string)

		//	var chunk_specifier_strings []string  /* elements like '23_51_3_1_135' */

		// read sequences from file into sequences slice
		// and set up maps from ids to indices and vice versa
		min_seq_len := 1000000000
		max_seq_len := -1

		index := 0
		scanner := bufio.NewScanner(fh)
		for scanner.Scan() {
			line := scanner.Text()
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
				sequences = append(sequences, line)
				fmt.Printf("   [%s]\n", line)
				index++
			}
		}
		if min_seq_len != max_seq_len {
			fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
			os.Exit(1)
		}
		seq_length := min_seq_len
		seq_chunk_set.Sequence_length = seq_length
		seq_chunk_set.Sequences = sequences
		seq_chunk_set.Index_id = index_id
		seq_chunk_set.Id_index = id_index
	} else {
		sequence_set = sequenceset.Construct_from_fasta_file(filename);
	}

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	var chunk_spec_strings []string
	var chunk_spec_arrays [][]int
	for used_so_far := 0; used_so_far+chunk_size <= seq_chunk_set.Sequence_length; used_so_far += chunk_size {
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
		fmt.Println(reflect.TypeOf(seq_matchindices))
		if !ok1 {
			seq_matchindices = make(map[string][]int)
			chunk__seq_matchindices[css] = seq_matchindices
		}
		for iseq, seq := range sequences { // loop over the sequences
			//	seq := sequences[iseq]
			//	seq_as_array := strings.Split(seq, "")
			fmt.Println("XXX: " + seq)
			chunk_seq := ""
			for j := 0; j < chunk_size; j += chunk_size { // assemble the sequence chunk
				chunk_seq += seq[csa[j] : csa[j]+1]
			}
			//	Chunk__seq_matchindices *map[string]*map[string]*[]int
			//	chunk__seq_matchindices := make(map[string]*map[string]*[]int)
			matchindices, ok2 := seq_matchindices[chunk_seq]
			if !ok2 {
				matchindices = make([]int, 1, 1)
				chunk__seq_matchindices[css][chunk_seq] = matchindices
			}
			chunk__seq_matchindices[css][chunk_seq] = append(matchindices, iseq) // push the sequence index onto appropriate slice
		}
	}

	seq_chunk_set.Chunk__seq_matchindices = chunk__seq_matchindices
	return &seq_chunk_set
}

func (scs sequence_chunk_set) Get_chunk_matchindex_counts(sequence string, matchindex_counts []int) []int {
	seq_length := len(sequence)
	if seq_length != scs.Sequence_length {
		fmt.Printf("sequence lengths: %8d %8d not equal; exiting.\n", seq_length, scs.Sequence_length)
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
					matchindex_counts[mindex]++
				}
			}
		}
	}
	n_seqs := len(matchindex_counts)
	sort.Ints(matchindex_counts)                 // small to large
	return matchindex_counts[n_seqs-10 : n_seqs] // sort return the last few ...
}
