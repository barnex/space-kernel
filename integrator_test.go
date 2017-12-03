package main

import (
	"math"
	"testing"
)

func TestIntegrator(t *testing.T) {
	cases := []struct {
		p0, v0       Vec
		dtheta       float64
		duration     float64
		pWant, vWant Vec
		tol          float64
	}{
		// tolerances drop quadratically with angle step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, Vec{0, -1}, 4e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, Vec{0, -1}, 4e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, Vec{0, -1}, 4e-7},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, Vec{0, -1}, 4e-9},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, Vec{0, -1}, 4e-11},
		// even smaller steps don't increase convergence anymore due to truncation.
	}

	for i, c := range cases {
		s := NewIntegrator(Gravity(1), c.dtheta)
		s.P = c.p0
		s.V = c.v0
		s.Advance(c.duration)
		if err := s.P.Sub(c.pWant).Len(); err > c.tol || math.IsNaN(err) {
			t.Errorf("case %v: p: have %v, want %v +/- %v, err=%v", i, s.P, c.pWant, c.tol, err)
		}
		if err := s.V.Sub(c.vWant).Len(); err > c.tol || math.IsNaN(err) {
			t.Errorf("case %v: v: have %v, want %v +/- %v, err=%v", i, s.V, c.vWant, c.tol, err)
		}
	}
}

func BenchmarkIntegrator(b *testing.B) {
	for i := 0; i < b.N; i++ {
		s := NewIntegrator(Gravity(1), 1e-3)
		s.P = Vec{1, 0}
		s.V = Vec{0, 1}
		s.Advance(2 * math.Pi)
	}
}

func BenchmarkIntegrate(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Integrate(AVerlet, Vec{1, 0}, Vec{0, 1}, Gravity(1), 0, 2*math.Pi, 1e-5, 1e-3)
	}
}

func TestAVerlet1(t *testing.T) {
	cases := []struct {
		p0, v0   Vec
		dtheta   float64
		duration float64
		pWant    Vec
		tol      float64
	}{
		// tolerances drop quadratically with angle step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 3e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 3e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 3e-7},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 3e-9},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 3e-11},
		// even smaller steps don't increase convergence anymore due to truncation.
	}

	stepper := AVerlet

	for i, c := range cases {
		p, _, _ := Integrate(stepper, c.p0, c.v0, Gravity(1), 0, c.duration, 1e-6, c.dtheta)
		if err := p.Sub(c.pWant).Len(); err > c.tol || math.IsNaN(err) {
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
		p, _, _ := Integrate(stepper, c.p0, c.v0, Gravity(1), 0, c.duration, c.dt, 0)
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
		// tolerances drop linearly with time step
		{Vec{0, 1}, Vec{1, 0}, 1e-1, math.Pi / 2, Vec{1, 0}, 2e-1},
		{Vec{0, 1}, Vec{1, 0}, 1e-2, math.Pi / 2, Vec{1, 0}, 2e-2},
		{Vec{0, 1}, Vec{1, 0}, 1e-3, math.Pi / 2, Vec{1, 0}, 2e-3},
		{Vec{0, 1}, Vec{1, 0}, 1e-4, math.Pi / 2, Vec{1, 0}, 2e-4},
		{Vec{0, 1}, Vec{1, 0}, 1e-5, math.Pi / 2, Vec{1, 0}, 2e-5},
		{Vec{0, 1}, Vec{1, 0}, 1e-6, math.Pi / 2, Vec{1, 0}, 2e-6},
	}

	stepper := SymEuler

	for i, c := range cases {
		p, _, _ := Integrate(stepper, c.p0, c.v0, Gravity(1), 0, c.duration, c.dt, 0)
		if err := p.Sub(c.pWant).Len(); err > c.tol || math.IsNaN(err) {
			t.Errorf("case %v: have %v, want %v +/- %v, err=%v", i, p, c.pWant, c.tol, err)
		}
	}
}
