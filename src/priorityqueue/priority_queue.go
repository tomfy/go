package priorityqueue

import (
//	"container/heap"
//	"fmt"
)

type IdCmf struct {
	Id   string
	Cmf float64
	Index int // The index of the item in the heap.
}

type PriorityQueue []*IdCmf

// implement the functions in heap interface
// this includes those of sort (Len, Less, and Swap)
func (pq PriorityQueue) Len() int { return len(pq) }

func (pq PriorityQueue) Less(i, j int) bool {
	// We want Pop to give us the item with the lowest Cmf, so use '<' here.
	return pq[i].Cmf < pq[j].Cmf
}

func (pq PriorityQueue) Swap(i, j int) {
	pq[i], pq[j] = pq[j], pq[i]
	pq[i].Index = i
	pq[j].Index = j
}

// and also Pop and Push

func (pq *PriorityQueue) Pop() interface{} {
	old := *pq
	n := len(old)
	item := old[n-1]
	item.Index = -1
	*pq = old[0 : n-1]
	return item
}

func (pq *PriorityQueue) Push(x interface{}) {
	n := len(*pq)
	item := x.(*IdCmf)
	item.Index = n
	*pq = append(*pq, item)
}

func (pq PriorityQueue) Peek() *IdCmf {
	return pq[0]
}

//



