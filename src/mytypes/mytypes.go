package mytypes

const MDchar = 'X' // the character used to indicate missing data, so e.g. 00X101X0211010X...

/*
type  Pair_int_int struct {
	A int
	B int
} /* */

type MatchInfo struct {
	Index              int
	Id                 string
	MatchCount         int
	OkChunkCount       int
	ChunkMatchFraction float64
}

type IdCmfDistance struct {
	Id                 string
	ChunkMatchFraction float64
	Distance           float64
}
