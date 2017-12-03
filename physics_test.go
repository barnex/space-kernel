package main

import (
	"testing"
)

// Test that the earth makes 1/2 orbit around the sun in 1/2 year.
// This verifies consistency of some constants,
// as well as the integrator working on this time/length scale.
func TestEarthOrbit(t *testing.T) {
	s := NewIntegrator(SunGravity(), 1e-3)
	start := Vec{EarthSMA, 0}
	s.P = start
	s.V = Vec{0, OrbitalV(SunMu, EarthSMA)}
	s.Advance(EarthP / 2)
	have := s.P
	want := Vec{-EarthSMA, 0}

	if have.Sub(want).Len() > EarthSMA/1e4 {
		t.Errorf("have: %v, want %v", have, want)
	}
}
