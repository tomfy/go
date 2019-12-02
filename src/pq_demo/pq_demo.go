package main

import (
	"container/heap"
	"fmt"
	"math/rand"
	"priorityqueue"
)

func main() {

	//	rand.Seed(seed)
	listItems := make([]*priorityqueue.IdCmf, 100, 100)
	//numbers = make([]float64, 0, 100)
	for i := 0; i < 100; i++ {
		number := float64(rand.Intn(100000)) * 0.00001
		x := priorityqueue.IdCmf{Id: "A_" + fmt.Sprintf("%5.5f", number), Cmf: number}
		listItems[i] = &x
	}


	// Initialize an empty priority queue
	priorityQueue := make(priorityqueue.PriorityQueue, 0)
	heap.Init(&priorityQueue)

	// Populate the pq
	for _, item := range listItems {
		pq_capped_push(&priorityQueue, item, 5)
	}

	// Print the order by Priority of expiry
	for priorityQueue.Len() > 0 {
	//	peek_item := (&priorityQueue).Peek()
		item := heap.Pop(&priorityQueue).(*priorityqueue.IdCmf)
	//	fmt.Printf("Id: %s Cmf: %5.3f \n", peek_item.Id, peek_item.Cmf)
		fmt.Printf("Id: %s Cmf: %5.5f \n", item.Id, item.Cmf)
	}

}


func pq_capped_push(pq *priorityqueue.PriorityQueue, x *priorityqueue.IdCmf, cap int) { // cap = 'capacity', max size.
	if len(*pq) < cap { // not at capacity, just add it
		heap.Push(pq, x)
	} else { // at capacity; compare to worst
		worst_cmf := pq.Peek().Cmf
		if x.Cmf > worst_cmf { // must bump one
			_ = heap.Pop(pq).(*priorityqueue.IdCmf) // discard worst one and
			heap.Push(pq, x) // add the new one.
		}
	}
}
