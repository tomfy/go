package mytypes

const MDchar = 'X' // the character used to indicate missing data, so e.g. 00X101X0211010X...

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
	Agmr               float64
	Hgmr               float64
	Distance           float64
}

type IdSeq struct {
	Id       string
	Sequence string
}

var Alpha float64 = 0.0
