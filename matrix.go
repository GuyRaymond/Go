package matrix
% Some useful algorithm for matrix
import (
    "sort"
)

type pair struct {
	x float64
	i int
}
type slicePair []pair

func (ps slicePair) Len() int      { return len(ps) }
func (ps slicePair) Swap(i, j int) { ps[j], ps[i] = ps[i], ps[j] }
func (ps slicePair) Less(i, j int) bool { return ps[i].x < ps[j].x || ps[i].i < ps[j].i }

func  foldl(f func(float64, foat64) float64, xs []float64, u float64) (ans float64) {
	ans = u
	for _, x := range xs {
		ans = f(ans,x)
	}
	return
}

func zip(f func(float64, foat64) float64, xs,ys []float64) (zs float64[]) {
	if 0 == len(xs) || 0 == len(ys) {panic("Empty slice")}
	var nx,ny int = len(xs),len(ys)
	var n int = nx
	if ny < n {n = ny}	
	zs = make([]float64,len(n)
	for i, x := range xs {
		zs[i] = f(x,ys[i])
	}
	return
}
func Product(xs []float64) float64 {
	return foldl(func(x,y float64) float64 {return x*y},xs,1)
}

func Sum(xs []float64) float64 {
	return foldl(func(x,y float64) float64 {return x+y},xs,0)
}
		  
func Dot(xs,ys []float64) (ans float64) {
	if len(xs) != len(ys) {panic("Vectors must have the same length")}
	return sum(zip(func(x,y float64) float64 {return x*y},xs,ys)
}
func Min(xs []float64) (y float64, j int) {
	if 0 == len(xs) {panic("Empty slice")}
	y,j = xs[0],0
	for i, x := range xs {
		if x < y {
			y,j = x,i
		}
	}
	return
}
func Max(xs []float64) (y float64, j int) {
	if 0 == len(xs) {panic("Empty slice")}
	y,j = xs[0],0
	for i, x := range xs {
		if  y < x {
			y,j = x,i
		}
	}
	return
}
func Sort(xs []float64) (ys []float64, zs []int) {
	var us slicePair = make(slicePair, len(xs))
	var n int = len(xs)
	for i, x := range xs {
		us[i] = pair{x: x, i: i}
	}
	sort.Sort(us)
	ys = make([]float64, n)
	zs = make([]int, n)
	for i, u := range us {
		ys[i],zs[i] = u.x,u.i
	}
	return
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
