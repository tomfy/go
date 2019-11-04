package sequenceset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	//	"sequence"
	"math/rand"
	"mytypes"
	"strings"
)

type Sequence_set struct {
	Sequences       []string
	Sequence_length int
	//	Seqs            []sequence.Sequence
	Index_id map[int]string
	Id_index map[string]int
}

func Construct_empty() *Sequence_set {
	var seq_set Sequence_set
	sequences := make([]string, 0)
	id_index := make(map[string]int)
	index_id := make(map[int]string)
	return &seq_set
}

func Construct_from_fasta_file(filename string) *Sequence_set {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	var seq_set Sequence_set

	var sequences []string
	//	var seqs []sequence.Sequence
	var id string
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
			id = match_strings[1]
			index_id[index] = id
			id_index[id] = index
			//		fmt.Println(match_strings[1])
		} else { // sequence line
			line = strings.TrimSpace(line)
			if len(line) < min_seq_len {
				min_seq_len = len(line)
			}
			if len(line) > max_seq_len {
				max_seq_len = len(line)
			}
			//	s := sequence.Construct(line, id)
			sequences = append(sequences, line)
			//	seqs = append(seqs, *s)
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
	//	seq_set.Seqs = seqs

	return &seq_set
}

func (set *Sequence_set) Add_sequence(index int, id string, sequence string) {
	if len(set.Sequences) > 0 {
		if len(sequence) != set.Sequence_length {
			os.Exit(1)
		}
	} else {
		set.Sequence_length = len(sequence)
	}

	set.Index_id[index] = id
	set.Id_index[id] = index
	set.Sequences = append(set.Sequences, sequence)
	//	fmt.Printf("%d %s  %d\n", index, id, len(set.Sequences))

}

func (set *Sequence_set) Add_missing_data(prob float64) {
	if prob > 0.0 {
		sequences := set.Sequences
		for i, _ := range sequences {
			chars := []byte(sequences[i])
			for j := 0; j < len(chars); j++ {
				if rand.Float64() < prob {
					// seq[i : i+1] = "X"
					chars[j] = mytypes.MDchar
				}
			}
			sequences[i] = string(chars)
			//	fmt.Println(sequences[i])
		}
		set.Sequences = sequences
	}
	return
}

func (set *Sequence_set) Index_to_id(index int) string {
	return set.Index_id[index]
}
