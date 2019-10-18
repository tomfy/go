package seqchunkset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"strings"
)

type sequence_chunk_set struct {
	Sequences *[]string
	Index_id  *map[int]string
	Id_index  *map[string]int
}

func Construct_from_fasta_file(filename string) *sequence_chunk_set {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	var seq_chunk_set sequence_chunk_set
	var sequences []string
	id_index := make(map[string]int)
	index_id := make(map[int]string)
	//	var chunk_specifier_strings []string  /* elements like '23_51_3_1_135' */
	min_seq_len := 1000000000
	max_seq_len := -1

	index := 0
	scanner1 := bufio.NewScanner(fh)
	for scanner1.Scan() {
		line := scanner1.Text()
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
	seq_chunk_set.Sequences = &sequences
	seq_chunk_set.Index_id = &index_id
	seq_chunk_set.Id_index = &id_index
	return &seq_chunk_set
}
