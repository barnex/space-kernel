package main

import (
	"math"
)

type Vec [2]float64

const (
	X = 0
	Y = 1
)

func (a Vec) Add(b Vec) Vec {
	return Vec{a[X] + b[X], a[Y] + b[Y]}
}

func (a Vec) MAdd(s float64, b Vec) Vec {
	return Vec{a[X] + s*b[X], a[Y] + s*b[Y]}
}

func (a Vec) Sub(b Vec) Vec {
	return Vec{a[X] - b[X], a[Y] - b[Y]}
}

func (a Vec) Dot(b Vec) float64 {
	return a[X]*b[X] + a[Y]*b[Y]
}

func (v Vec) Mul(a float64) Vec {
	return Vec{a * v[X], a * v[Y]}
}

func (v Vec) Mul3(a Vec) Vec {
	return Vec{a[X] * v[X], a[Y] * v[Y]}
}

func (v Vec) Div(a float64) Vec {
	return v.Mul(1 / a)
}

// Length (norm).
func (v Vec) Len() float64 {
	return math.Sqrt(v[X]*v[X] + v[Y]*v[Y])
}

// Length squared
func (v Vec) Len2() float64 {
	return v[X]*v[X] + v[Y]*v[Y]
}

// Returns a copy of v, scaled to unit length.
func (v Vec) Normalized() Vec {
	l := 1 / math.Sqrt(v[X]*v[X]+v[Y]*v[Y])
	return Vec{v[X] * l, v[Y] * l}
}
