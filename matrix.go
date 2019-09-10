package matrix

// Some useful algorithm for matrix
import (
	"math"
	"sort"
)
//A = [...]float64{{0.8147, 0.0975, 0.1576}, {0.9058, 0.2785, 0.9706}, {0.1270, 0.5469, 0.9572}, 
// {0.9134, 0.9575, 0.4854}, {0.6324, 0.9649, 0.8003}}
type pair struct {
	x float64
	i int
}
type slicePair []pair

func (ps slicePair) Len() int           { return len(ps) }
func (ps slicePair) Swap(i, j int)      { ps[j], ps[i] = ps[i], ps[j] }
func (ps slicePair) Less(i, j int) bool { return ps[i].x < ps[j].x || ps[i].i < ps[j].i }


func Fill(xs []float64, x float64) {
	for i := 0; i < len(xs); i++ { xs[i] = x }
}

func FillMatrix(xs [][]float64, x float64) {
	for i := 0; i < len(xs); i++ { Fill(xs[i],x) }
}

func Eye(int m,n) (xs [][]float64) {
	xs = make([][]float64,m)
	for i := 0; i < m; i++ {
		xs[i] = make([]float64,n)
		xs[i][i] = 1
	}
	return
}

func Ones(int m,n) (xs [][]float64) {
	xs = make([][]float64,m)
	for i := 0; i < m; i++ {
		xs[i] = make([]float64,n)
		for j := 0; i < n; j++ { xs[i][j] = x }	
	}
	return
}

func Size(xs [][]float64) (rows, cols int) {
	rows = len(xs)
	if 0 == rows {
		cols = 0
	} else {
		cols = len(xs[0])
		var c int
		for i := 1; i < rows; i++ {
			c = len(xs[i])
			if c < cols {
				cols = c
			}
		}
	}
	return
}

func Col(xs [][]float64, j float64) (ys []float64) {
	var m,_ int = Size(xs)
	ys = make([]float64,m)
	for i := 0; i < m; i++ { ys[i] = xs[i][j] }
}

func ColMatrix(xs [][]float64, j float64) (ys [][]float64) {
	var m,n int = Size(xs)
	ys = make([][]float64,m)
	for i := 0; i < m; i++ { 
		ys[i] = = make([]float64,1)
		ys[i][0] = xs[i][j] 
	}
}

func Row(xs [][]float64, j float64) (ys []float64) {
	var _,n int = Size(xs)
	ys = make([]float64,n)
	for j := 0; j < n; j++ { ys[i] = xs[i][j] }
}

func rowMatrix(xs [][]float64, j float64) (ys [][]float64) {
	var m,n int = Size(xs)
	ys = make([][]float64,n)
	for j := 0; j < m; j++ { 
		ys[j] = make([]float64,1)
		ys[j][0] = xs[i][j] 
	}
}

func Transpose(xs [][]float64) (ys [][]float64) {
	rows, cols = Size(xs)
	ys = make([][]float64, rows)
	if 0 < rows {
		cols = len(xs[0])
		var c int
		for i := 0; i < rows; i++ {
			ys[i] = make([][]float64, cols)
			for j := 0; j < cols; j++ {
				ys[i][j] = xs[j][i]
			}
		}
	}
	return
}

func Map(f func(float64) float64, xs []float64) (ys []float64) {
	var n int = len(xs)
	ys = make([]float64, n)
	for i, x := range xs {
		ys[i] = f(x)
	}
	return
}

func MapArray(f func([]float64) []float64, xs [][]float64) (ys [][]float64) {
	var n int = len(xs)
	ys = make([]float64, n)
	for i, x := range xs {
		ys[i] = f(x)
	}
	return
}

func Apply(f func(float64) float64, xs [][]float64) (ys [][]float64) {
	var n int = len(xs)
	ys = make([][]float64, n)
	for i, _ := range xs {
		ys[i] = Map(f,xs[i])
	}
	return
}

func Filter(f func(float64) bool, xs []float64) (ys []float64) {
	var n int = len(xs)
	ys = make([]float64, n)
	var j int
	for i, x := range xs {
		if f(x) {
			ys[j] = x
			j++
		}
	}
	ys = ys[0:j]
	return
}

func FilterArray(f func([]float64) bool, xs [][]float64) (ys [][]float64) {
	var n int = len(xs)
	ys = make([]float64, n)
	var j int
	for i, x := range xs {
		if f(x) {
			ys[j] = make([]float64, len(x))
			copy(ys[j],x)
			j++
		}
	}
	ys = ys[0:j]
	return
}

func foldl(f func(float64, foat64) float64, xs []float64, u float64) (ans float64) {
	ans = u
	for _, x := range xs {
		ans = f(ans, x)
	}
	return
}

func foldlMatrix(f func([]float64, []foat64) []float64, xs [][]float64, u []float64) (ans []float64) {
	ans = u
	for _, x := range xs {
		ans = f(ans, x)
	}
	return
}

func zip(f func(float64, foat64) float64, xs, ys []float64) (zs []float64) {
	if 0 == len(xs) || 0 == len(ys) {
		panic("Empty slice")
	}
	var nx, ny int = len(xs), len(ys)
	var n int = nx
	if ny < n {
		n = ny
	}
	zs = make([]float64, n)
	for i, x := range xs {
		zs[i] = f(x, ys[i])
	}
	return
}

func zipMatrix(f func([]float64, []foat64) []float64, xs, ys [][]float64) (zs [][]float64) {
	if 0 == len(xs) || 0 == len(ys) {
		panic("Empty slice")
	}
	var nx, ny int = len(xs), len(ys)
	var n int = nx
	if ny < n {
		n = ny
	}
	zs = make([][]float64, n)
	for i := o; i < n; i++ {
		zs[i] = zip(f, xs[i], ys[i])
	}
	return
}

func Plus(xs, ys []float64) []float64 {
	if len(xs) != len(ys) {
		panic("Vectors must have the same length")
	}
	return zip(func(x, y float64) float64 { return x + y }, xs, ys)
}

func PlusMatrix(xs, ys [][]float64) (zs [][]float64) {
	if len(xs) != len(ys) {
		panic("Vectors must have the same number of rows")
	}
	return zipMatrix(Plus, xs, ys)
}

func Minus(xs, ys []float64) []float64 {
	if len(xs) != len(ys) {
		panic("Vectors must have the same length")
	}
	return zip(func(x, y float64) float64 { return x - y }, xs, ys)
}

func MinusMatrix(xs, ys [][]float64) (zs [][]float64) {
	if len(xs) != len(ys) {
		panic("Vectors must have the same number of rows")
	}
	return zipMatrix(Minus, xs, ys)
}

func Product(xs []float64) float64 {
	return foldl(func(x, y float64) float64 { return x * y }, xs, 1)
}

func Sum(xs []float64) float64 {
	return foldl(func(x, y float64) float64 { return x + y }, xs, 0)
}

func SumCol(xs [][]float64) float64 {
	return foldl(Plus, xs, 0)
}

func SumRow(xs [][]float64) float64 {
	return SumCol(Transpose(xs))
}

func SumArray(xs [][]float64) float64 {
	return Sum(SumCol(xs))
}

func NormSquared(xs []float64) (s float64) {
	s = 0
	for _, x := range xs {
		s += x*x
	}
	return
}

func Norm(xs []float64) {
	return math.Sqrt(NormSquared(xs))
}

func NormArray(xs [][]float64) (s float64) {
	s = 0
	for _, x := range xs {
		s += NormSquared(x)
	}
	s = math.Sqrt(s)
	return
}

func Minus(xs, ys []float64) []float64 {
	if len(xs) != len(ys) {
		panic("Vectors must have the same length")
	}
	return zip(func(x, y float64) float64 { return x - y }, xs, ys)
}

func Min(xs ...float64) (y float64, j int) {
	if 0 == len(xs) {
		panic("Empty slice")
	}
	y, j = xs[0], 0
	for i, x := range xs {
		if x < y {
			y, j = x, i
		}
	}
	return
}

func Max(xs ...float64) (y float64, j int) {
	if 0 == len(xs) {
		panic("Empty slice")
	}
	y, j = xs[0], 0
	for i, x := range xs {
		if y < x {
			y, j = x, i
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
		ys[i], zs[i] = u.x, u.i
	}
	return
}

func Mult(xs, ys [][]float64) (zs [][]float64) {
	mx, nx := Size(xs)
	my, ny := Size(ys)
	if nx != my {
		panic("number of columns of the left matrix must equal to the number of rows of the right matrix")
	}
	zs = make([][]float64, mx)
	for i,x := range xs {
		zs[i] = Sumcol(Apply(func(u float64) float64 { return x * u}, ys)
	}
	return
}

func DotArray(xs, ys [][]float64) ([][]float64) {
	return Mult(Transpose(xs),ys)
}


func householderqr(A [][]float64) (Q,R [][]float64)  {
// QR de decomposition with the method of Householder reflections
	m,n := size(A)
	Q = Eye(m,m)
	R = make([][]float64,m)
	for i := 0; i < m; i++ {
		R[i] = make([]float64,n)
		copy(R[i],A[i])
	}
	var a,rho,nu float64
	var u []float64
	var us [][1]float64
	for k := 0; k < n; k++ {
		u = R[k:m][k]
		us = make([][1]float64,m-k+1)
		rho = -1.0
		if u[0] < 0 {rho = 1}
		nu = Norm(u)
		a = rho*Norm(u)
		u[0] = u[0] - a
		u = Map(func(x float64) float64 { return x/nu },u)
		v := Mult(Map(func(x float64) float64 { return 2*x },u),transpose(u))
		R(k:m,k:n) = R(k:m,k:n) - v*R(k:m,k:n)
		if  1 == k  {
			Q(k:m,k:m) = Q(k:m,k:m) - v*Q(k:m,k:m)
		} else  {
			p = k-1
			Q(k:m,1:p) = Q(k:m,1:p) - v*Q(k:m,1:p)
			Q(k:m,k:m) = Q(k:m,k:m) - v*Q(k:m,k:m)
		}
	}
	Q = transpose(Q)
	return
}
