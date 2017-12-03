package main

import "math"

// Some constants related to our solar system and universe.
const (
	// gravitational constant in m3/kgs2
	// https://en.wikipedia.org/wiki/Gravitational_constant
	G = 6.674e-11

	// G * Sun mass
	// en.wikipedia.org/wiki/Standard_gravitational_parameter
	SunMu = 1.32712440042e20

	// https://en.wikipedia.org/wiki/Moon
	MoonSMA  = 384399000       // semi-major axis in m
	MoonP    = 27.321661 * Day // orbital period in s
	MoonSynP = 29.530589 * Day // synodic period in s

	// en.wikipedia.org/wiki/Earth
	EarthMu  = 3.986004418e14      // gravitational parameter
	EarthSMA = 149598023000        // semi-major axis in m
	EarthP   = 365.256363004 * Day //  orbital period in s

	Day = 24 * 3600 // earth day in seconds
)

// Gravity returns the gravitational acceleration
// around a point mass with gravitational parameter mu (M*G).
func Gravity(mu float64) AccelFunc {
	return func(p Vec, t float64) Vec {
		x, y := p[X], p[Y]
		len2 := x*x + y*y
		len3 := len2 * math.Sqrt(len2)
		f := -mu / len3
		return Vec{f * x, f * y}
	}
}

// SunGravity returns the gravitational acceleration of the sun,
// located at the origin.
func SunGravity() AccelFunc {
	return Gravity(SunMu)
}

// EarthGravity returns the gravitational acceleration of the earth,
// located at the origin.
func EarthGravity() AccelFunc {
	return Gravity(EarthMu)
}

// OrbitalV returns the orbital velocity of an object
// orbiting with semi-major axis sma around an object with gravitational parameter mu.
func OrbitalV(mu, sma float64) float64 {
	return math.Sqrt(mu / sma)
}
