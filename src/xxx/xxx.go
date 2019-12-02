package main

import (
	"fmt"
	"sync"
	"time"
)

var wg1, wg2 sync.WaitGroup

type Id_d struct {
	Id       string
	Distance float64
}

func pinger(c chan Id_d, a Id_d) {
	defer wg2.Done()
		c <- a
}

func printer(c chan string) {
	for {
		msg := <-c
		fmt.Println(msg)
		time.Sleep(time.Second * 1)
	}
}

func load_into_map(c chan Id_d, m *map[string]float64) {
	defer wg1.Done()
	for {
		idd, ok := <-c
		if !ok {
			return
		}
		(*m)[idd.Id] = idd.Distance
	}
}

func main() {

	var ch chan Id_d = make(chan Id_d)
	the_map := make(map[string]float64)
	wg1.Add(1) // add 1 to the count of goroutines whose completion we are waiting for.
	go load_into_map(ch, &the_map)
	for i := 0; i < 30; i++ {
		id := "A_" + fmt.Sprintf("%d", i)
		d := float64(i) * 0.01234
		x := Id_d{id, d}

		wg2.Add(1)
		go pinger(ch, x)
	}
	wg2.Wait()  // wait until pinger goroutines have all returned.
	close(ch)
	wg1.Wait() // wait until load_into_map goroutine returns.

	fmt.Println(len(the_map))
	for id, d := range the_map {
		fmt.Printf("%s  %8.4f\n", id, d)
	}
}
