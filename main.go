package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"os"
)

var (
	flagR     = flag.Float64("r", 1, "starting distance from center")
	flagV     = flag.Float64("v", 1, "starting velocity")
	flagT     = flag.Float64("t", 356*Day, "total time")
	flagDTh   = flag.Float64("d", 1e-3, "angular resolution")
	flagEvery = flag.Float64("e", 1, "output every `N` steps")
)

func main() {
	flag.Parse()
	log.SetPrefix("#")

	Mtotal := 100.0
	MFuel := 90.0
	Ve := 3000.0           // m/s
	Exhaust := 1.0         // kg/s
	Thrust := Ve * Exhaust // kgm/s2 == N

	accel := func(_ Vec, t float64) Vec {
		Mex := Exhaust * t
		if Mex > MFuel {
			return Vec{}
		}
		M := Mtotal - Mex
		return Vec{0, Thrust / M}
	}

	s := NewIntegrator(accel, 1e-4)
	s.OutputEvery = 1
	s.Advance(100)

	Mfinal := Mtotal - MFuel
	want := Ve * math.Log(Mtotal/Mfinal)
	log.Println(want)
}

func Log(x ...interface{}) {
	fmt.Fprintln(os.Stderr, x...)
}

func min(a, b float64) float64 {
	if a < b {
		return a
	}
	return b
}
