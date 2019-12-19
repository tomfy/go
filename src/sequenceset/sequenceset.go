package sequenceset

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	//	"sequence"
	//	"math"
	"math/rand"
	"mytypes"
	"priorityqueue"
	"sort"
	"strconv"
	"strings"
	"sync"
	"etc"
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

type QsetSsetQSmi struct {
	Qss   *Sequence_set
	Sss   *Sequence_set
	Qs_mi [][]*mytypes.MatchInfo
}

func Construct_from_matrix_file(filename string, max_md_prop float64, id_seqset *map[string]*Sequence_set,
	seq_set *Sequence_set, waitgroup *sync.WaitGroup) {

	defer waitgroup.Done()
	
	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
//	var seq_set Sequence_set
	var sequences []string

	id_index := make(map[string]int)
	index_id := make(map[int]string)
	marker_id_index := make(map[string]int)
	marker_index_id := make(map[int]string)

	// read sequences from file into sequences slice
	// and set up maps from ids to indices and vice versa
	min_seq_len := 1000000000
	max_seq_len := -1

	seq_index := 0
	scanner := bufio.NewScanner(fh)
	scanner.Buffer(make([]byte, 10000), 1000000) // th
	line_number := 0
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line) // split on one or more whitespace chars.
		//	fmt.Fprintln(os.Stderr, "line number: ", line_number)
		if line_number == 0 { // this line should have MARKER and then marker ids (tab separated)
			if fields[0] != "MARKER" {
				os.Exit(1)
			}
			marker_ids := fields[1:] // everything after "MARKER"
			for index, id := range marker_ids {
				marker_id_index[id] = index
				marker_index_id[index] = id
			}

		} else {
			seq_id := fields[0]
			index_id[seq_index] = seq_id
			id_index[seq_id] = seq_index
			genotypes_sequence := strings.Join(fields[1:], "")
			//	fmt.Fprintln(os.Stderr, "g sequence: ", genotypes_sequence)
			sequences = append(sequences, genotypes_sequence)
			if len(genotypes_sequence) < min_seq_len {
				min_seq_len = len(genotypes_sequence)
			}
			if len(genotypes_sequence) > max_seq_len {
				max_seq_len = len(genotypes_sequence)
			}
			seq_index++
		} /* */
		line_number++
	} // end of reading line lines of file
	if min_seq_len != max_seq_len { // all sequence lengths must be the same, otherwise exit.
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}
	seq_length := min_seq_len
	seq_set.Sequence_length = seq_length
	seq_set.Sequences = sequences
	seq_set.SeqIndex_id = index_id
	seq_set.SeqId_index = id_index

	for id, _ := range id_index {
		(*id_seqset)[id] = seq_set
	}

	seq_set.SnpIndex_id = marker_index_id // make(map[int]string)
	seq_set.SnpId_index = marker_id_index // make(map[string]int)
	//	fmt.Fprintln(os.Stderr, "before missing_data_counts")
	seq_set.missing_data_counts()
	//	fmt.Fprintln(os.Stderr, "after missing_data_counts")
	max_md_count := int(max_md_prop * float64(seq_length))
	seq_set.Max_md_count = max_md_count

	//   sort by amount of missing data
	seq_set.Sorted_snp_ids, seq_set.N_ok_snps = keys_sorted_by_value(seq_set.SnpId_mdcount, seq_set.Max_md_count)
	/*	for i, snp_id := range seq_set.Sorted_snp_ids {
			fmt.Println(i, snp_id, seq_set.SnpId_mdcount[snp_id])
	} /* */

//	return &seq_set
}

func Construct_sets_from_matrix_file(filename string, n_sets_to_make int, max_md_prop float64,
	id_seqset *map[string]*Sequence_set, seq_sets []*Sequence_set) {

	fh, err := os.Open(filename)
	if err != nil {
		os.Exit(1)
	}
//	var seq_set Sequence_set
	var sequences []string

	id_index := make(map[string]int)
	index_id := make(map[int]string)
	marker_id_index := make(map[string]int)
	marker_index_id := make(map[int]string)

	// read sequences from file into sequences slice
	// and set up maps from ids to indices and vice versa
	min_seq_len := 1000000000
	max_seq_len := -1

	seq_index := 0
	scanner := bufio.NewScanner(fh)
	scanner.Buffer(make([]byte, 10000), 1000000) // th
	line_number := 0
	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line) // split on one or more whitespace chars.
		//	fmt.Fprintln(os.Stderr, "line number: ", line_number)
		if line_number == 0 { // this line should have MARKER and then marker ids (tab separated)
			if fields[0] != "MARKER" {
				os.Exit(1)
			}
			marker_ids := fields[1:] // everything after "MARKER"
			for index, id := range marker_ids {
				marker_id_index[id] = index
				marker_index_id[index] = id
			}

		} else {
			seq_id := fields[0]
			index_id[seq_index] = seq_id
			id_index[seq_id] = seq_index
			genotypes_sequence := strings.Join(fields[1:], "")
			//	fmt.Fprintln(os.Stderr, "g sequence: ", genotypes_sequence)
			sequences = append(sequences, genotypes_sequence)
			if len(genotypes_sequence) < min_seq_len {
				min_seq_len = len(genotypes_sequence)
			}
			if len(genotypes_sequence) > max_seq_len {
				max_seq_len = len(genotypes_sequence)
			}
			seq_index++
		} /* */
		line_number++
	} // end of reading lines of file
	if min_seq_len != max_seq_len { // all sequence lengths must be the same, otherwise exit.
		fmt.Printf("min, max sequence lengths: %8d %8d lengths not equal; exiting.\n", min_seq_len, max_seq_len)
		os.Exit(1)
	}
	seq_length := min_seq_len

	n_sequences := len(sequences)

	n_seqs_in_each_set := n_sequences/n_sets_to_make + 1

//	seq_sets := make([]*Sequence_set, n_sets_to_make)
	n_seqs_used_so_far := 0
	for i := 0; i < n_sets_to_make; i++ {
		a_seq_set := Sequence_set{}
		set_size := etc.MinInt(n_seqs_in_each_set, n_sequences-n_seqs_used_so_far)

		set_seq_id_index := make(map[string]int)
		set_seq_index_id := make(map[int]string)
		set_marker_id_index := make(map[string]int)
		set_marker_index_id := make(map[int]string)

		first_seq_index := n_seqs_used_so_far
		a_seq_set.Sequences = sequences[first_seq_index : first_seq_index+set_size] //
		for j := 0; j < set_size; j++ {                                             // get the seq_id_index, etc. for each set
			overall_seq_index := n_seqs_used_so_far + j
			seq_id := index_id[overall_seq_index]
			set_seq_id_index[seq_id] = j
			set_seq_index_id[j] = seq_id

			(*id_seqset)[seq_id] = &a_seq_set
		}
		for marker_id, marker_index := range marker_id_index { // get the marker_id_index, etc. for each set
			set_marker_id_index[marker_id] = marker_index
			set_marker_index_id[marker_index] = marker_id
		}
		n_seqs_used_so_far += set_size

		a_seq_set.SeqId_index = set_seq_id_index
		a_seq_set.SeqIndex_id = set_seq_index_id
		a_seq_set.SnpId_index = set_marker_id_index
		a_seq_set.SnpIndex_id = set_marker_index_id
		a_seq_set.Sequence_length = seq_length

		a_seq_set.missing_data_counts()
		//	fmt.Fprintln(os.Stderr, "after missing_data_counts")
		max_md_count := int(max_md_prop * float64(seq_length))
		a_seq_set.Max_md_count = max_md_count

		//   sort by amount of missing data
		a_seq_set.Sorted_snp_ids, a_seq_set.N_ok_snps = keys_sorted_by_value(a_seq_set.SnpId_mdcount, a_seq_set.Max_md_count)

		seq_sets[i] = &a_seq_set
	}

//	return seq_sets
}

func Construct_from_fasta_file(filename string, max_md_prop float64, rand_md_rate float64, id_seqset *map[string]*Sequence_set) *Sequence_set {
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
		r, _ := regexp.Compile(`^>(\S+)`)
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

	for id, _ := range id_index {
		(*id_seqset)[id] = &seq_set
	}

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
		fmt.Println(i, snp_id, seq_set.SnpId_mdcount[snp_id])
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

func (set *Sequence_set) missing_data_counts() { // for the snp with id snp_id, count the number of sequences with missing data
	//	fmt.Fprintf(os.Stderr, "nnn snps: %d\n", len(set.SnpId_index))
	set.SnpId_mdcount = make(map[string]int)
//	fmt.Fprintln(os.Stderr, "len SnpId_index: ", len(set.SnpId_index))
	for snp_id, snp_idx := range set.SnpId_index {
		set.SnpId_mdcount[snp_id] = 0
		for _, seq := range set.Sequences { // loop over sequences
			//	fmt.Println("BBB: ", seq)
			if seq[snp_idx:snp_idx+1] == string(mytypes.MDchar) {
				//	fmt.Println("snp_id: ", snp_id, snp_idx, seq[snp_idx:snp_idx+1])
				set.SnpId_mdcount[snp_id]++
			}

		}
		//	fmt.Fprintln(os.Stderr, "snp id, md count: ", snp_id, set.SnpId_mdcount[snp_id])
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

// /*
func (seq_set *Sequence_set) Fraction_of_all_distances_AA(prob float64) map[string]map[string]float64 {

	idpair_dist := make(map[string]map[string]float64)
	for i1, seq1 := range seq_set.Sequences {
		id1 := seq_set.Seq_index_to_id(i1)

		for i2 := i1 + 1; i2 < len(seq_set.Sequences); i2++ {
			if prob >= 1 || rand.Float64() < prob {
				id2 := seq_set.Seq_index_to_id(i2)
				seq2 := seq_set.Sequences[i2]
				//	dist_old := distance_old(sseq, qseq)
				n00_22, n11, nd1, nd2 := Distance(seq1, seq2)
				dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
				//	fmt.Printf("%v  %v\n", dist_old, dist)
				id2_dist, ok := idpair_dist[id1]
				if !ok {
					id2_dist = make(map[string]float64)
				}
				id2_dist[id2] = dist
				idpair_dist[id1] = id2_dist

				id1_dist, ok := idpair_dist[id2]
				if !ok {
					id1_dist = make(map[string]float64)
				}
				id1_dist[id1] = dist
				idpair_dist[id2] = id1_dist

			}

		}
	}
	return idpair_dist
} /* */

func (seq_set *Sequence_set) Fraction_of_all_distances_AB(seq_set2 *Sequence_set, prob float64) map[string]map[string]float64 {

	idpair_dist := make(map[string]map[string]float64)
	for i1, seq1 := range seq_set.Sequences {
		id1 := seq_set.Seq_index_to_id(i1)

		for i2, seq2 := range seq_set2.Sequences {
			if prob >= 1 || rand.Float64() < prob {
				id2 := seq_set2.Seq_index_to_id(i2)
				//	seq2 := seq_set2.Sequences[i2]
				//	dist_old := distance_old(sseq, qseq)
				n00_22, n11, nd1, nd2 := Distance(seq1, seq2)
				dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
				//	fmt.Printf("%v  %v\n", dist_old, dist)
				id2_dist, ok := idpair_dist[id1]
				if !ok {
					id2_dist = make(map[string]float64)
					id2_dist[id2] = dist
					idpair_dist[id1] = id2_dist
				} else {
					id2_dist[id2] = dist
				}

				/*	id1_dist, ok := idpair_dist[id2]
					if !ok {
						id1_dist = make(map[string]float64)
					}
					id1_dist[id1] = dist
					idpair_dist[id2] = id1_dist /* */

			}

		}
	}
	return idpair_dist
} /* */

func (q_seq_set *Sequence_set) Candidate_distances_qs(s_seq_set *Sequence_set, qid_matchcandidates map[string][]*mytypes.MatchInfo, qid_matches map[string][]mytypes.IdCmfDistance) int {
	// for the candidate matches in qid_matchcandidates, get the full distances, sort, and output.
	// this is for the case of searching for matches between two distinct Sequence_sets (i.e. 'AB')
	dist_calc_count := 0
	for qindex, qseq := range q_seq_set.Sequences {
		qid := q_seq_set.Seq_index_to_id(qindex)
		matchcandidates := qid_matchcandidates[qid]

		id_matchcount_distance_triples := make([]mytypes.IdCmfDistance, len(matchcandidates))
		fmt.Printf("%s   ", qid)
		for i, matchinfo := range matchcandidates {
			sseq_index := matchinfo.Index
			sseq_id := matchinfo.Id // s_seq_set.Seq_index_to_id(sseq_index)
			sseq := s_seq_set.Sequences[sseq_index]
			//	dist_old := distance_old(sseq, qseq)
			n00_22, n11, nd1, nd2 := Distance(sseq, qseq)
			dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
			//	hgmr := float64(nd2) / float64(n00_22 + nd2)
			//	agmr := float64(nd1+nd2) / float64(n00_22+n11+nd1+nd2)
			//	fmt.Printf("%v  %v\n", dist_old, dist)

			dist_calc_count++
			matchinfo := mytypes.IdCmfDistance{sseq_id, matchinfo.ChunkMatchFraction, dist}
			id_matchcount_distance_triples[i] = matchinfo
			qid_matches[qid] = append(qid_matches[qid], matchinfo)
		}

		sort.Slice(id_matchcount_distance_triples,
			func(i, j int) bool {
				return id_matchcount_distance_triples[i].Distance < id_matchcount_distance_triples[j].Distance
			})
		// output the best matches to stdout:1
		for _, a_triple := range id_matchcount_distance_triples {
			fmt.Printf("%s %6.5f %6.5f  ", a_triple.Id, a_triple.ChunkMatchFraction, a_triple.Distance)
		}
		fmt.Printf("\n")
	}
	return dist_calc_count
} /* */

func (q_seq_set *Sequence_set) Candidate_distances_pq(s_seq_set *Sequence_set, qid_cmfpq map[string]*priorityqueue.PriorityQueue) /*, qid_matches map[string][]mytypes.IdCmfDistance)*/ int {
	// for the candidate matches in qid_matchcandidates, get the full distances, sort, and output.
	// this is for the case of searching for matches between two distinct Sequence_sets (i.e. 'AB')
	dist_calc_count := 0
	for qindex, qseq := range q_seq_set.Sequences {
		qid := q_seq_set.Seq_index_to_id(qindex)
		matchcandidates := qid_cmfpq[qid]

		id_matchcount_distance_triples := make([]mytypes.IdCmfDistance, len(*matchcandidates))
		fmt.Printf("%s   ", qid)
		for i, matchinfo := range *matchcandidates {
			// sseq_index := matchinfo.Index

			sseq_id := matchinfo.Id // s_seq_set.Seq_index_to_id(sseq_index)
			sseq_index, ok := s_seq_set.SeqId_index[sseq_id]
			if !ok {
				os.Exit(1)
			}
			sseq := s_seq_set.Sequences[sseq_index]
			//	dist_old := distance_old(sseq, qseq)
			n00_22, n11, nd1, nd2 := Distance(sseq, qseq)
			dist := float64(nd1+2*nd2) / float64(n00_22+n11+nd1+nd2)
			//	hgmr := float64(nd2) / float64(n00_22 + nd2)
			//	agmr := float64(nd1+nd2) / float64(n00_22+n11+nd1+nd2)
			//	fmt.Printf("%v  %v\n", dist_old, dist)
			/*if(dist != distx){
				os.Exit(1)
			}*/
			dist_calc_count++
			matchinfo := mytypes.IdCmfDistance{sseq_id, matchinfo.Cmf, dist} // ChunkMatchFraction, dist}
			id_matchcount_distance_triples[i] = matchinfo
			//	qid_matches[qid] = append(qid_matches[qid], matchinfo)
		}

		sort.Slice(id_matchcount_distance_triples,
			func(i, j int) bool {
				return id_matchcount_distance_triples[i].Distance < id_matchcount_distance_triples[j].Distance
			})

		for _, a_triple := range id_matchcount_distance_triples {
			fmt.Printf("%s %6.5f %6.5f  ", a_triple.Id, a_triple.ChunkMatchFraction, a_triple.Distance)
		}
		fmt.Printf("\n")
	}
	return dist_calc_count
}

// *************** not methods ***********************

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
//	fmt.Fprintf(os.Stderr, "max_mds %d  n_ok: %d \n", max_mds, n_ok)
	return keys, n_ok // the keys sorted by value (small to large), and the number of values <= max_mds
}

func Distance(seq1 string, seq2 string) (int, int, int, int) {
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

/* func min(a int, b int) int {
	if a < b {
		return a
	}
	return b
} /* */

func (seq_set *Sequence_set) Check_seq_index_id_maps() bool {
	for i, _ := range seq_set.Sequences {
		id := seq_set.SeqIndex_id[i]
		ii := seq_set.SeqId_index[id]
		if ii != i {
			fmt.Println("seq set index/id inconsistency. i, ii: ", i, ii)
			os.Exit(1)
		}
	}
	for id, idx := range seq_set.SeqId_index {
		id2 := seq_set.SeqIndex_id[idx]
		if id2 != id {
	fmt.Println("seq set index/id inconsistency. id, id2: ", id, id2)
			os.Exit(1)
		}
	}
	return true
}
