package qseqchunkset

import (
	"container/heap"
	"fmt"
	"math/rand"
	"mytypes"
	"os"
	"sort"
	"sync"

	"priorityqueue"
	"sequenceset"
)

// var Match_count_increment_count int

type Q_sequence_chunk_set struct {
	Sequence_set              *sequenceset.Sequence_set
	Chunk_size                int
	Chunk_specs               []chunk_spec
	Missing_data_chunk_counts []int // Missing_data_chunk_counts[i] is the number of chunks with missing data for sequence with index i
	N_chunked_sequences       int   // number of sequences entered into Chunk__seq_matchindices, Missing_data_chunk_counts so far.

	Chunk__seq_matchindices   map[string]map[string][]int // keys are chunk specifier strings, values are maps, whose keys are chunk seqs (e.g. '10020101') values are slices of sequence indices
	Id__chsp_chseq            map[string]map[string]string
}

type chunk_spec struct {
	s string   // e.g. 'A0 A34 A5 A101'
	a []string //
}

/* type Index_matchcount struct {
	index      int
	matchcount int
} /* */

// *********************** exported funct
