package main

import "math"

func NewIntegrator(a AccelFunc, dtheta float64) *Integrator {
	return &Integrator{
		acc:    a,
		dtheta: dtheta,
		dt:     1e-5, // TODO
	}
}

type Integrator struct {
	T    float64
	P, V Vec

	acc    AccelFunc
	dtheta float64
	dt     float64
}

func (s *Integrator) Advance(duration float64) {
	p, v := s.P, s.V
	dt := s.dt
	dth := s.dtheta
	t := s.T
	tmax := t + duration

	a1 := s.acc(p, t)
	var dt2 float64
	for t+dt < tmax {

		dt_2 := dt / 2
		dt2_2 := dt * dt_2

		p[X] += dt*v[X] + dt2_2*a1[X]
		p[Y] += dt*v[Y] + dt2_2*a1[Y]

		a2 := s.acc(p, t)
		v[X] += dt_2 * (a1[X] + a2[X])
		v[Y] += dt_2 * (a1[Y] + a2[Y])

		// do a second step
		p[X] += dt*v[X] + dt2_2*a2[X]
		p[Y] += dt*v[Y] + dt2_2*a2[Y]
		a1 = s.acc(p, t+dt)
		v[X] += dt_2 * (a1[X] + a2[X])
		v[Y] += dt_2 * (a1[Y] + a2[Y])

		//da := math.Sqrt(a1.Sub(a2).Len2() / (0.25 * a1.Add(a2).Len2()))
		da := math.Sqrt((sqr(a1[X]-a2[X]) + sqr(a1[Y]-a2[Y])) / (0.25 * (sqr(a1[X]+a2[X]) + sqr(a1[Y]+a2[Y]))))
		fac := clamp(dth / da)
		dt2 = dt * fac

		//a1 = a2 // FSAL
		t += dt * 2
		dt = dt2
	}
	s.dt = dt2

	dt = tmax - t
	if dt != 0 {
		p, v, _ = AVerlet(p, v, s.acc, t, dt, dth)
		t += dt
	}

	s.P = p
	s.V = v
	s.T = t
}
func sqr(x float64) float64 { return x * x }

// Integrate advances position and velocity given an acceleration function,
// stepping from t0 to tmax with initial time step dt0 and angle step d theta.
func Integrate(s StepperFunc, p0, v0 Vec, a AccelFunc, t0, tmax, dt0, dth float64) (Vec, Vec, float64) {
	t := t0
	dt := dt0
	p := p0
	v := v0

	var dt2 float64
	for t+dt < tmax {
		p, v, dt2 = s(p, v, a, t, dt, dth)
		t += dt
		dt = dt2
	}
	dt = tmax - t
	if dt != 0 {
		p, v, _ = s(p, v, a, t, dt, dth)
		t += dt
	}

	return p, v, dt2
}

// An AccelFunc function returns acceleration (m/s2) as a function of position and time.
type AccelFunc func(p Vec, t float64) Vec

// A StepperFunc function makes one integration step
// and returns the new position, velocity and time step.
type StepperFunc func(p, v Vec, acc AccelFunc, t, dt, maxDa float64) (Vec, Vec, float64)

// Adaptive Verlet stepper.
func AVerlet(p, v Vec, acc AccelFunc, t, dt, maxDa float64) (Vec, Vec, float64) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2

	a1 := acc(p, t)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := acc(p, t)
	v = v.MAdd(dt_2, a1.Add(a2))

	da := a1.Sub(a2).Len() / (0.5 * a1.Add(a2).Len())

	fac := clamp(maxDa / da)
	dt *= fac

	return p, v, dt
}

// Velocity Verlet stepper.
// https://en.wikipedia.org/wiki/Verlet_integration
func Verlet(p, v Vec, acc AccelFunc, t, dt, _ float64) (Vec, Vec, float64) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2
	a1 := acc(p, t)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := acc(p, t)
	v = v.MAdd(dt_2, a1.Add(a2))
	return p, v, dt
}

// Symplectic Euler stepper, used for debugging.
// https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
func SymEuler(p, v Vec, acc AccelFunc, t, dt, _ float64) (Vec, Vec, float64) {
	v = v.MAdd(dt, acc(p, t))
	p = p.MAdd(dt, v)
	return p, v, dt
}

func clamp(x float64) float64 {
	if x < 0.5 {
		return 0.5
	}
	if x > 2 {
		return 2
	}
	return x
}
