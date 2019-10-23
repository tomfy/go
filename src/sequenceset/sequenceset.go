package sequenceset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"strings"
)

type Sequence_set struct {
	Sequences               []string
	Sequence_length         int
	Index_id                map[int]string
	Id_index                map[string]int
}

func Construct_from_fasta_file(filename string) *Sequence_set {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	var seq_set Sequence_set

	var sequences []string
	id_index := make(map[string]int)
	index_id := make(map[int]string)

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
	//		fmt.Println(match_strings[1])
		} else { // sequence line
			line = strings.TrimSpace(line)
			if len(line) < min_seq_len {
				min_seq_len = len(line)
			}
			if len(line) > max_seq_len {
				max_seq_len = len(line)
			}
			sequences = append(sequences, line)
	//		fmt.Printf("   [%s]\n", line)
			index++
		}
	}
	if min_seq_len != max_seq_len { // all sequence lengths must be the same, otherwise exit.
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}
	seq_length := min_seq_len
	seq_set.Sequence_length = seq_length
	seq_set.Sequences = sequences
	seq_set.Index_id = index_id
	seq_set.Id_index = id_index

	return &seq_set
}

func (set Sequence_set) Index_to_id(index int) string {
	return set.Index_id[index]
}
