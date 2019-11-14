package sequenceset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	//	"sequence"
	"sort"
	"math/rand"
	"mytypes"
	"strings"
	"strconv"
)

type Sequence_set struct {
	Sequences       []string
	Sequence_length int
	SeqIndex_id     map[int]string
	SeqId_index     map[string]int
	SnpIndex_id     map[int]string
	SnpId_index     map[string]int
	SnpId_mdcount   map[string]int // keys: snp ids, values: missing data counts.
	Sorted_snp_ids  []string // sorted by missing data, low to high.
	Max_md_count int // only use snps with <= this number of seqs with missing data
	N_ok_snps int // the number of snps with sufficiently small amount of missing data
	
}

func Construct_empty() *Sequence_set {
	seq_set := Sequence_set{
		make([]string, 0),    // Sequences
		-1,                   // Sequence_length
		make(map[int]string), // SeqIndex_id
		make(map[string]int), // SeqId_index
		make(map[int]string), // SnpIndex_id
		make(map[string]int), // SnpId_index
		make(map[string]int), // SnpId_mdcount
		make([]string, 0),    // Sorted_snp_ids
		-1,                   // Max_md_count
		-1,                   // N_ok_snps
		
	}
	return &seq_set
}

func Construct_from_fasta_file(filename string, max_md_prop float64, rand_md_rate float64) *Sequence_set {
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
	var seq_set Sequence_set
	var sequences []string
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
		r, _ := regexp.Compile("^>([A-Za-z0-9]+).*")
		if r.MatchString(line) { // >id line
			match_strings := r.FindStringSubmatch(line)
			id = match_strings[1]
			index_id[index] = id
			id_index[id] = index
			// fmt.Println(match_strings[1])
		} else { // sequence line
			line = strings.TrimSpace(line)
			if len(line) < min_seq_len {
				min_seq_len = len(line)
			}
			if len(line) > max_seq_len {
				max_seq_len = len(line)
			}
			// s := sequence.Construct(line, id)
			sequences = append(sequences, line)
			// seqs = append(seqs, *s)
			// fmt.Printf("   [%s]\n", line)
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
	seq_set.SeqIndex_id = index_id
	seq_set.SeqId_index = id_index

	seq_set.SnpIndex_id = make(map[int]string)
	seq_set.SnpId_index = make(map[string]int)
	
	seq_set.add_snp_ids() //	
	seq_set.Add_missing_data(rand_md_rate)


	seq_set.SnpId_mdcount = make(map[string]int)
	for _, snp_id := range seq_set.SnpIndex_id {
		// seq_set.SnpId_mdcount[snp_id] =
			seq_set.missing_data_count(snp_id)	
	}

/*	for snpid, mdcount := range seq_set.SnpId_mdcount{
		fmt.Println(snpid, mdcount)
	} */
	
	// sort by amount of missing data	
	seq_set.Sorted_snp_ids = keys_sorted_by_value(seq_set.SnpId_mdcount)
	for i, snp_id := range seq_set.Sorted_snp_ids {
		fmt.Println(i, snp_id, seq_set.SnpId_mdcount[snp_id])
	}
	
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

	set.SeqIndex_id[index] = id
	set.SeqId_index[id] = index
	set.Sequences = append(set.Sequences, sequence)
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
		//	fmt.Println(string(chars))
			sequences[i] = string(chars)
			//	fmt.Println("AAAA", sequences[i])
		}
		set.Sequences = sequences
	}
	return
}

func (set *Sequence_set) Seq_index_to_id(index int) string {
	return set.SeqIndex_id[index]
}

/* func (set *Sequence_set) Prune_snps(max_missing_data_proportion) {

	pruned_seqs := make([]string, len(set.Sequences))
	pruned_seq_length_so_far := 0
		for j, g := chars {
			md_prop := float64(set.Missing_data_seq_counts[j]) / len(set.Sequences)
			if md_prop <= max_missing_data_proportion {
				for i, seq := range set.Sequences {
					pruned_seqs[i] += seq[
					chars := []byte(seq)

				pruned_chars = append(pruned_chars, g)

}
		}
	}
	return set.SeqIndex_id[index]
} /* */

func (set *Sequence_set) missing_data_count(snp_id string) { // for the snp with id snp_id, count the number of sequences with missing data
	snp_idx := set.SnpId_index[snp_id]
	for _, seq := range set.Sequences { // loop over sequences
	//	fmt.Println("BBB: ", seq)
		if seq[snp_idx:snp_idx+1] == string(mytypes.MDchar) {
			fmt.Println("snp_id: ", snp_id, snp_idx, seq[snp_idx: snp_idx+1])
			set.SnpId_mdcount[snp_id]++
		}
	}
}

// /*
func (set *Sequence_set) add_snp_ids() { // for now, just make an id by prepending 'A' to the position within sequence as it appears in fasta file.
	for i := 0; i < set.Sequence_length; i++ {
		snp_id := "A" + strconv.Itoa(i)
		set.SnpIndex_id[i] = snp_id
		set.SnpId_index[snp_id] = i
	}
} /* */

// *************** not methods ***********************

/* func count_missing_data_snps_in_seq(seq string) int {
	count := 0
	for j := 0; j < len(seq); j++ { // loop over snps
		if seq[j] == mytypes.MDchar {
			count++
		}
	}
	return count
} /* 

func count_seqs_w_missing_snp_data(seq string) int {
	count := 0
	for j := 0; j < len(seq); j++ { // loop over snps
		if seq[j] == mytypes.MDchar {
			count++
		}
	}
	return count
} /**/

// for a map w string keys, int values, return a slice containing the keys, but sorted with smallest value ones at beginning.
func keys_sorted_by_value(amap map[string]int) []string {
	keys := make([]string, 0, len(amap))
	for k := range amap {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return amap[keys[i]] < amap[keys[j]] } )
	return keys
}
