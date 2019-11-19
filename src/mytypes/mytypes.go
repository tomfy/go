package mytypes

const MDchar = 'X' // the character used to indicate missing data, so e.g. 00X101X0211010X...

type  Pair_int_int struct {
	A int
	B int
}

type IntIntIntF64 struct {
	A int
	B int
	C int
	D float64
}

type StringF64F64 struct {
	A string
	B float64
	C float64
}
