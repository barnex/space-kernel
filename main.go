package main

import (
	"flag"
	"fmt"
)

var flagVerlet = flag.Bool("v", false, "use verlet")

func main() {
	flag.Parse()

	solver := SymEuler
	if *flagVerlet {
		solver = Verlet
	}

	//every := 0.01
	max := 3.

	for h := 0.5; h > 5e-8; h /= 2 {
		p := Vec{1, 0, 0}
		v := Vec{0, 1.2, 0}
		//for t := 0.0; t < max; t += every {
		//	p, v = solver(p, v, every, h)
		//	fmt.Println(t, p[X], p[Y], h)
		//}
		p, v = solver(p, v, max, h)
		fmt.Println(h, p[X], p[Y], v[X], v[Y])
	}
}

func Acc(p Vec) Vec {
	return p.Normalized().Div(-p.Len2())
}

// https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
func SymEuler(p, v Vec, d, dt float64) (Vec, Vec) {
	for t := 0.0; t < d; t += dt {
		v = v.MAdd(dt, Acc(p))
		p = p.MAdd(dt, v)
	}
	return p, v
}

// https://en.wikipedia.org/wiki/Verlet_integration
func Verlet(p, v Vec, d, dt float64) (Vec, Vec) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2
	for t := 0.0; t < d; t += dt {
		a1 := Acc(p)
		p = p.MAdd(dt, v).MAdd(dt2_2, a1)
		a2 := Acc(p)
		v = v.MAdd(dt_2, a1.Add(a2))
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
