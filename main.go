package main

import (
	"flag"
	"fmt"
)

func main() {
	flag.Parse()

	solver := AVerlet

	//every := 0.01
	max := 900.
	h := 1e-3

	p := Vec{1, 0, 0}
	v := Vec{0, 1.4, 0}
	//i := 0
	for t := 0.0; t < max; t += h {
		p, v, h = solver(p, v, h)
		//i++
		//if i == 5 {
		fmt.Println(t, p[X], p[Y], h)
		//i = 0
		//}
	}
}

func Acc(p Vec) Vec {
	return p.Normalized().Div(-p.Len2())
}

// https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
func SymEuler(p, v Vec, dt float64) (Vec, Vec, float64) {
	v = v.MAdd(dt, Acc(p))
	p = p.MAdd(dt, v)
	return p, v, dt
}

// https://en.wikipedia.org/wiki/Verlet_integration
func Verlet(p, v Vec, dt float64) (Vec, Vec, float64) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2
	a1 := Acc(p)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := Acc(p)
	v = v.MAdd(dt_2, a1.Add(a2))
	return p, v, dt
}

// Adaptive verlet
func AVerlet(p, v Vec, dt float64) (Vec, Vec, float64) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2

	a1 := Acc(p)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := Acc(p)
	v = v.MAdd(dt_2, a1.Add(a2))

	da := a1.Sub(a2).Len() / (a1.Add(a2).Len())
	maxDa := 1e-2

	fac := clamp(maxDa / da)
	dt *= fac

	return p, v, dt
}

func clamp(x float64) float64 {
	if x < 0.1 {
		return 0.1
	}
	if x > 10 {
		return 10
	}
	return x
}

//type Rocket struct {
//	Pos     Vec     // position in m
//	V       Vec     // velocity in m/s
//	Heading Vec     // nose direction (unit vector)
//	M       float64 // mass in kg
//	Thrust  float64 // Thrust in N
//	Vsp     float64 // specific velocity (specific impulse*g)
//
