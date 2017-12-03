package main

import (
	"flag"
	"fmt"
	"math"
	"os"
)

var (
	flagR     = flag.Float64("r", 1, "starting distance from center")
	flagV     = flag.Float64("v", 1, "starting velocity")
	flagT     = flag.Float64("t", Earth.P, "total time")
	flagDTh   = flag.Float64("d", 1e-3, "angular resolution")
	flagEvery = flag.Float64("e", 1, "output every `N` steps")
)

func main() {
	flag.Parse()

	max := *flagT
	h := 1.0
	dth := *flagDTh
	every := *flagEvery

	p := Vec{Earth.R + MoonSMA, 0}
	MoonV := 2 * math.Pi * MoonSMA / MoonP
	v := Vec{0, Earth.V() + MoonV}

	i := every
	n := 0
	for t := 0.0; t < max; {
		p, v, h = AVerlet(p, v, AccSolar, t, h, dth)
		t += h
		if i == every {
			e := Earth.Pos(t)
			m := MoonPos(t)
			fmt.Println(t, p[X], p[Y], m[X], m[Y], e[X], e[Y], h)
			i = 0
		}
		i++
		n++
	}

	Log(n, "steps")
	Log("dt", h, "s")
}

type Planet struct {
	Mu float64 // M*G in m3/s2
	R  float64 // Semi-major axis in m
	P  float64 // Orbital period in s
}

var (
	Earth = Planet{
		Mu: EarthMu,
		R:  EarthSMA,
		P:  EarthP,
	}
)

func (b *Planet) V() float64 {
	return 2 * math.Pi * b.R / b.P
}

func (p *Planet) Pos(t float64) Vec {
	return Vec{
		X: p.R * math.Cos(2*math.Pi*t/p.P),
		Y: p.R * math.Sin(2*math.Pi*t/p.P),
	}
}

func MoonPos(t float64) Vec {
	m := Vec{
		X: MoonSMA * math.Cos(2*math.Pi*t/MoonP),
		Y: MoonSMA * math.Sin(2*math.Pi*t/MoonP),
	}
	return Earth.Pos(t).Add(m)
}

func (p *Planet) Acc(r Vec, t float64) Vec {
	l := r.Sub(p.Pos(t))
	return l.Normalized().Mul(-p.Mu / l.Len2())
}

func UnitAcc(p Vec, t float64) Vec {
	return p.Normalized().Mul(-1 / p.Len2())
}

func AccSolar(p Vec, t float64) Vec {
	return AccSun(p).Add(Earth.Acc(p, t))
}

func AccSun(p Vec) Vec {
	return p.Normalized().Mul(-SunMu / p.Len2())
}

func Log(x ...interface{}) {
	fmt.Fprintln(os.Stderr, x...)
}

//type Rocket struct {
//	Pos     Vec     // position in m
//	V       Vec     // velocity in m/s
//	Heading Vec     // nose direction (unit vector)
//	M       float64 // mass in kg
//	Thrust  float64 // Thrust in N
//	Vsp     float64 // specific velocity (specific impulse*g)
//
