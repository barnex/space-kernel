package main

import (
	"math"
)

type Vec [3]float64

const (
	X = 0
	Y = 1
	Z = 2
)

func (a Vec) Add(b Vec) Vec {
	return Vec{a[X] + b[X], a[Y] + b[Y], a[Z] + b[Z]}
}

func (a Vec) MAdd(s float64, b Vec) Vec {
	return Vec{a[X] + s*b[X], a[Y] + s*b[Y], a[Z] + s*b[Z]}
}

func (a Vec) Sub(b Vec) Vec {
	return Vec{a[X] - b[X], a[Y] - b[Y], a[Z] - b[Z]}
}

func (a Vec) Dot(b Vec) float64 {
	return a[X]*b[X] + a[Y]*b[Y] + a[Z]*b[Z]
}

func (v Vec) Mul(a float64) Vec {
	return Vec{a * v[X], a * v[Y], a * v[Z]}
}

func (v Vec) Mul3(a Vec) Vec {
	return Vec{a[X] * v[X], a[Y] * v[Y], a[Z] * v[Z]}
}

func (v Vec) Div(a float64) Vec {
	return v.Mul(1 / a)
}

// Length (norm).
func (v Vec) Len() float64 {
	return math.Sqrt(v.Dot(v))
}

// Length squared
func (v Vec) Len2() float64 {
	return v.Dot(v)
}

// Returns a copy of v, scaled to unit length.
func (v Vec) Normalized() Vec {
	l := 1 / math.Sqrt(v[X]*v[X]+v[Y]*v[Y]+v[Z]*v[Z])
	return Vec{v[X] * l, v[Y] * l, v[Z] * l}
}

func (a Vec) Cross(b Vec) Vec {
	x := a[Y]*b[Z] - a[Z]*b[Y]
	y := a[Z]*b[X] - a[X]*b[Z]
	z := a[X]*b[Y] - a[Y]*b[X]
	return Vec{x, y, z}
}
