package main

import (
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
		// tolerances drop quadratically with time step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 3e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 3e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 3e-7},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 3e-9},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 3e-11},
		// even smaller steps don't increase convergence anymore due to truncation.
	}

	stepper := AVerlet

	for i, c := range cases {
		p, _ := Integrate(stepper, c.p0, c.v0, UnitAcc, 0, c.duration, 1e-6, c.dtheta)
		if err := p.Sub(c.pWant).Len(); err > c.tol {
			t.Errorf("case %v: have %v, want %v +/- %v, err=%v", i, p, c.pWant, c.tol, err)
		}
	}
}

func TestVerlet1(t *testing.T) {
	cases := []struct {
		p0, v0   Vec
		dt       float64
		duration float64
		pWant    Vec
		tol      float64
	}{
		// tolerances drop quadratically with time step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 3e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 3e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 3e-7},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 3e-9},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 3e-11},
		// even smaller steps don't increase convergence anymore due to truncation.
	}

	stepper := Verlet

	for i, c := range cases {
		p, _ := Integrate(stepper, c.p0, c.v0, UnitAcc, 0, c.duration, c.dt, 0)
		if err := p.Sub(c.pWant).Len(); err > c.tol {
			t.Errorf("case %v: have %v, want %v +/- %v, err=%v", i, p, c.pWant, c.tol, err)
		}
	}
}

func TestSymEuler1(t *testing.T) {
	cases := []struct {
		p0, v0   Vec
		dt       float64
		duration float64
		pWant    Vec
		tol      float64
	}{
		// tolerances drop linerarly with time step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 2e-1},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 2e-2},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 2e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 2e-4},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 2e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-6, math.Pi / 2, Vec{1, 0}, 2e-6},
		// even smaller steps don't increase convergence anymore due to truncation.
	}

	stepper := SymEuler

	for i, c := range cases {
		p, _ := Integrate(stepper, c.p0, c.v0, UnitAcc, 0, c.duration, c.dt, 0)
		if err := p.Sub(c.pWant).Len(); err > c.tol {
			t.Errorf("case %v: have %v, want %v +/- %v, err=%v", i, p, c.pWant, c.tol, err)
		}
	}
}
