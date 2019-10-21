package seqchunkset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"strings"
)

type sequence_chunk_set struct {
	Sequences               *[]string
	Sequence_length         int
	Chunk_size              int
	Chunk_spec_strings      *[]string
	Chunk_spec_arrays       *[]*[]int
	Index_id                *map[int]string
	Id_index                *map[string]int
	Chunk__seq_matchindices *map[string]*map[string]int
}

func Construct_from_fasta_file(filename string, chunk_size int) *sequence_chunk_set {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	var seq_chunk_set sequence_chunk_set
	var sequences []string
	id_index := make(map[string]int)
	index_id := make(map[int]string)
	chunk__seq_matchindices := make(map[string]*map[string]int)
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
	seq_chunk_set.Sequence_length = min_seq_len
	seq_chunk_set.Sequences = &sequences
	seq_chunk_set.Index_id = &index_id
	seq_chunk_set.Id_index = &id_index

	// get the chunk-specifying strings (e.g. "1_5_109_4") and store in a slice
	var chunk_spec_strings []string
	var chunk_spec_arrays []*[]int
	for used_so_far := 0; used_so_far+chunk_size <= seq_chunk_set.Sequence_length; used_so_far += chunk_size {
		chunk_spec_str := fmt.Sprintf("%d", used_so_far)
		var chunk_spec_array []int
		chunk_spec_array = append(chunk_spec_array, used_so_far)
		for i := 1; i < chunk_size; i++ {
			chunk_spec_str += fmt.Sprintf("_%d", used_so_far+i)
			chunk_spec_array = append(chunk_spec_array, used_so_far+i)
		}
		chunk_spec_strings = append(chunk_spec_strings, chunk_spec_str)
		chunk_spec_arrays = append(chunk_spec_arrays, &chunk_spec_array)
	}
	seq_chunk_set.Chunk_spec_strings = &chunk_spec_strings
	seq_chunk_set.Chunk_spec_arrays = &chunk_spec_arrays

	// for each sequence, consider as chunks and store sequence index in
	for iseq := range sequences {
		seq := sequences[iseq]
		for ich := range chunk_spec_arrays {
			csa := chunk_spec_arrays[ich] // a pointer to a slice
			
		}
	}

	return &seq_chunk_set
}
