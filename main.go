package main

import "fmt"

func main() {

	p := Vec{1, 0, 0}
	v := Vec{0, 1.2, 0}

	h := 1e-3
	t := 0.0
	for t < 1000 {
		e := 0.0
		for e < .1 {
			v = v.MAdd(h, Acc(p))
			p = p.MAdd(h, v)
			e += h
		}
		t += e
		fmt.Println(t, p[X], p[Y])
	}
}

func Acc(p Vec) Vec {
	return p.Normalized().Div(-p.Len2())
}

func Euler(p, v Vec, d, dt float64) (Vec, Vec) {
	for t := 0.0; t < d; t += dt {
		v = v.MAdd(h, Acc(p))
		p = p.MAdd(h, v)
	}
	return p, v
}

//type Rocket struct {
//	Pos     Vec     // position in m
//	V       Vec     // velocity in m/s
//	Heading Vec     // nose direction (unit vector)
//	M       float64 // mass in kg
//	Thrust  float64 // Thrust in N
//	Vsp     float64 // specific velocity (specific impulse*g)
//
