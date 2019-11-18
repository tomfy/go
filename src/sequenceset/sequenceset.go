package sequenceset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	//	"sequence"
	"math/rand"
	"mytypes"
	"sort"
	"strconv"
	"strings"
)

type Sequence_set struct {
	Sequences       []string
	Sequence_length int
	SeqIndex_id     map[int]string
	SeqId_index     map[string]int
	SnpIndex_id     map[int]string
	SnpId_index     map[string]int
	SnpId_mdcount   map[string]int // keys: snp ids, values: missing data counts.
	Sorted_snp_ids  []string       // sorted by missing data, low to high.
	Max_md_count    int            // only use snps with <= this number of seqs with missing data
	N_ok_snps       int            // the number of snps with sufficiently small amount of missing data

}

/* func Construct_empty() *Sequence_set {
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
} /* */

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

	seq_set.add_snp_ids() // for now these are just A0, A1, A2, etc
	fmt.Fprintf(os.Stderr, "n snps: %d %d\n", len(seq_set.SnpId_index), len(seq_set.SnpIndex_id))

	seq_set.Add_missing_data(rand_md_rate)
	fmt.Fprintf(os.Stderr, "n snps: %d %d\n", len(seq_set.SnpId_index), len(seq_set.SnpIndex_id))


//	for snp_index, snp_id := range seq_set.SnpIndex_id {
		// seq_set.SnpId_mdcount[snp_id] =
//		fmt.Fprintf(os.Stderr, "snp index, id: %d  %s\n", snp_index, snp_id)
		seq_set.missing_data_counts()
//	}

	max_md_count := int(max_md_prop * float64(seq_length))
	seq_set.Max_md_count = max_md_count
	/*	for snpid, mdcount := range seq_set.SnpId_mdcount{
		fmt.Println(snpid, mdcount)
	} */

	// sort by amount of missing data
	seq_set.Sorted_snp_ids, seq_set.N_ok_snps = keys_sorted_by_value(seq_set.SnpId_mdcount, seq_set.Max_md_count)
 /*	for i, snp_id := range seq_set.Sorted_snp_ids {
	//	fmt.Println(i, snp_id, seq_set.SnpId_mdcount[snp_id])

	} /* */

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

func (set *Sequence_set) missing_data_counts() { // for the snp with id snp_id, count the number of sequences with missing data
	//	fmt.Fprintf(os.Stderr, "nnn snps: %d\n", len(set.SnpId_index))
	set.SnpId_mdcount = make(map[string]int)
	for snp_id, snp_idx := range set.SnpId_index {
		set.SnpId_mdcount[snp_id] = 0
		for _, seq := range set.Sequences { // loop over sequences
			//	fmt.Println("BBB: ", seq)
			if seq[snp_idx:snp_idx+1] == string(mytypes.MDchar) {
			//	fmt.Println("snp_id: ", snp_id, snp_idx, seq[snp_idx:snp_idx+1])
				set.SnpId_mdcount[snp_id]++
			}
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

func (seq_set *Sequence_set) Candidate_distances_AA(qid_matchcandidates  map[string][]*mytypes.IntIntIntF64) {
 // for the candidate matches in qid_matchcandidates, get the full distances, sort, and output.
 // this is for the case of searching for matches within a single data Sequence_set (i.e. 'AA') 
	for qindex, qseq := range seq_set.Sequences {
			qid := seq_set.Seq_index_to_id(qindex)
			matchcandidates := qid_matchcandidates[qid]

			id_matchcount_distance_triples := make([]mytypes.StringF64F64, len(matchcandidates))
			fmt.Printf("%s   ", qid)
			for i, matchinfo := range matchcandidates {
				sseq_index := matchinfo.A
				sseq_id := seq_set.Seq_index_to_id(sseq_index)
				sseq := seq_set.Sequences[sseq_index]
			//	dist_old := distance_old(sseq, qseq)
				n00_22, n11, nd1, nd2 := distance(sseq, qseq)
				dist := float64(nd1 + 2*nd2)/float64(n00_22 + n11 + nd1 + nd2)
			//	fmt.Printf("%v  %v\n", dist_old, dist)
				/*if(dist != distx){
					os.Exit(1)
				}*/
				id_matchcount_distance_triples[i] = mytypes.StringF64F64{sseq_id, matchinfo.D, dist}
			}
			sort.Slice(id_matchcount_distance_triples,
				func(i, j int) bool { return id_matchcount_distance_triples[i].C < id_matchcount_distance_triples[j].C })

			for _, a_triple := range id_matchcount_distance_triples {
				fmt.Printf("%s %6.5f %6.5f  ", a_triple.A, a_triple.B, a_triple.C)
			}
			fmt.Printf("\n")
		}

}

func (q_seq_set *Sequence_set) Candidate_distances_AB(s_seq_set *Sequence_set, qid_matchcandidates  map[string][]*mytypes.IntIntIntF64) {
 // for the candidate matches in qid_matchcandidates, get the full distances, sort, and output.
 // this is for the case of searching for matches within a single data Sequence_set (i.e. 'AA') 
	for qindex, qseq := range q_seq_set.Sequences {
			qid := q_seq_set.Seq_index_to_id(qindex)
			matchcandidates := qid_matchcandidates[qid]

			id_matchcount_distance_triples := make([]mytypes.StringF64F64, len(matchcandidates))
			fmt.Printf("%s   ", qid)
			for i, matchinfo := range matchcandidates {
				sseq_index := matchinfo.A
				sseq_id := s_seq_set.Seq_index_to_id(sseq_index)
				sseq := s_seq_set.Sequences[sseq_index]
			//	dist_old := distance_old(sseq, qseq)
				n00_22, n11, nd1, nd2 := distance(sseq, qseq)
				dist := float64(nd1 + 2*nd2)/float64(n00_22 + n11 + nd1 + nd2)
			//	fmt.Printf("%v  %v\n", dist_old, dist)
				/*if(dist != distx){
					os.Exit(1)
				}*/
				id_matchcount_distance_triples[i] = mytypes.StringF64F64{sseq_id, matchinfo.D, dist}
			}
			sort.Slice(id_matchcount_distance_triples,
				func(i, j int) bool { return id_matchcount_distance_triples[i].C < id_matchcount_distance_triples[j].C })

			for _, a_triple := range id_matchcount_distance_triples {
				fmt.Printf("%s %6.5f %6.5f  ", a_triple.A, a_triple.B, a_triple.C)
			}
			fmt.Printf("\n")
		}

}

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
func keys_sorted_by_value(amap map[string]int, max_mds int) ([]string, int) {
	keys := make([]string, 0, len(amap))
	n_ok := 0
	for k, v := range amap {
		keys = append(keys, k)
	//	fmt.Fprintf(os.Stderr, "k v: %s  %d\n", k, v)
		if v <= max_mds {
			n_ok++
		}
	}
	sort.Slice(keys, func(i, j int) bool { return amap[keys[i]] < amap[keys[j]] })
	fmt.Fprintf(os.Stderr, "max_mds %d  n_ok: %d \n", max_mds, n_ok)
	return keys, n_ok // the keys sorted by value (small to large), and the number of values <= max_mds
}


func distance(seq1 string, seq2 string) (int, int, int, int) {
	//	zero_count := 0
	one_count := 0
	two_count := 0
	n00_22 := 0 // homozygous, no change, i.e. 0->0 or 2->2
	n11 := 0    // heterozygous, no change, i.e. 1 -> 1
	// n02 := 0 // same as two_count
	// n01_12 := 0 // same as one_count
	for i := 0; i < len(seq1); i++ {
		c1 := seq1[i : i+1]
		c2 := seq2[i : i+1]
		if c1 == "0" {
			if c2 == "0" {
				//	zero_count++
				n00_22++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				two_count++
			}
		} else if c1 == "1" {
			if c2 == "0" {
				one_count++
			} else if c2 == "1" {
				//	zero_count++
				n11++
			} else if c2 == "2" {
				one_count++
			}
		} else if c1 == "2" {
			if c2 == "0" {
				two_count++
			} else if c2 == "1" {
				one_count++
			} else if c2 == "2" {
				//	zero_count++
				n00_22++
			}
		}
	}
	return n00_22, n11, one_count, two_count
	/*	ok_count := n00_22 + n11 + one_count + two_count // number of sites where neither seq has missing data
		dist_count := one_count + 2*two_count          // sums differences, i.e. 0-1 -> +=1, 0-2 -> += 2, ...
		var distance float64
		if ok_count > 0 {
			distance = float64(dist_count) / float64(ok_count)
		} else {
			distance = -1.0 // couldn't calculate because no sites without missing data
		}
		return distance */
}
