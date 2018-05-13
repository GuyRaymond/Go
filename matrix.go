package matrix
% Some useful algorithm for matrix
import (
    gosort"sort"
)

type pair struct {
	x float64
	i int
}
type slicePair []pair

func (ps slicePair) Len() int      { return len(ps) }
func (ps slicePair) Swap(i, j int) { ps[j], ps[i] = ps[i], ps[j] }
func (ps slicePair) Less(i, j int) bool { return ps[i].x < ps[j].x || ps[i].i < ps[j].i }

func sum(xs []float64) float64 {
    val ans float64 = 0;
    for _,x := range xs {
        ans += x
    }
    return ans
}

func product(xs []float64) float64 {
    val ans float64 = 1;
    for _,x := range xs {
        ans *= x
    }
    return ans
}


func dot(x,y []float64) float64 {
}

func sort(xs []float64) (ys []float64, zs []int) {
	var us slicePair = make(slicePair, len(xs))
	var n int = len(xs)
	for i, x := range xs {
		us[i] = pair{x: x, i: i}
	}
	gosort.Sort(us)
	ys = make([]float64, n)
	zs = make([]int, n)
	for i, u := range us {
		ys[i],zs[i] = u.x,u.i
	}
	return
}
