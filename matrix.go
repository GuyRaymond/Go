package matrix
% Some useful algorithm for matrix
import (
    "sort"
)

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

func sort(xs []float64) (float64[],uint[]) {
    ys := make [][2]float64
    for i,x := range xs {
    }
}
