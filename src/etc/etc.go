package etc

// various useful functions ...

func MinInt(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func MaxInt(a int, b int) int {
        if a > b {
                return a
        }
        return b
}

