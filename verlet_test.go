package main

import (
	"log"
	"math"
	"testing"
)

func TestAVerlet1(t *testing.T) {
	cases := []struct {
		p0, v0   Vec
		dtheta   float64
		duration float64
		pWant    Vec
		tol      float64
	}{
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 0},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 0},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 0},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 0},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 0},
	}

	for i, c := range cases {
		p, _ := Integrate(c.p0, c.v0, UnitAcc, 0, c.duration, 1e-6, c.dtheta)
		log.Println(p)
		log.Println(c.pWant)
		log.Println(p.Sub(c.pWant).Len())
		if err := p.Sub(c.pWant).Len(); err > c.tol {
			t.Errorf("case %v: have %v, want %v +/- %v, err=%v", i, p, c.pWant, c.tol, err)
		}
	}
}
