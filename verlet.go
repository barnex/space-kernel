package main

// Adaptive Verlet integrator.
func AVerlet(p, v Vec, acc Accel, t, dt, maxDa float64) (Vec, Vec, float64) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2

	a1 := acc(p, t)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := acc(p, t)
	v = v.MAdd(dt_2, a1.Add(a2))

	da := a1.Sub(a2).Len() / (a1.Add(a2).Len())

	fac := clamp(maxDa / da)
	dt *= fac

	return p, v, dt
}

func Integrate(p0, v0 Vec, a Accel, t0, tmax, dt0, dth float64) (Vec, Vec) {

	t := t0
	dt := dt0
	p := p0
	v := v0

	for t+dt < tmax {
		p, v, dt = AVerlet(p, v, a, t, dt, dth)
		//p, v = Verlet(p, v, a, t, dt)
		t += dt
	}
	// shorter final step to reach exact tmax
	dt = tmax - t
	p, v, _ = AVerlet(p, v, a, t, dt, dth)
	t += dt

	return p, v
}

type Accel func(p Vec, t float64) Vec

func clamp(x float64) float64 {
	if x < 0.5 {
		return 0.5
	}
	if x > 2 {
		return 2
	}
	return x
}

// Velocity Verlet integrator.
// https://en.wikipedia.org/wiki/Verlet_integration
func Verlet(p, v Vec, acc Accel, t, dt float64) (Vec, Vec) {
	dt2_2 := dt * dt / 2
	dt_2 := dt / 2
	a1 := acc(p, t)
	p = p.MAdd(dt, v).MAdd(dt2_2, a1)
	a2 := acc(p, t)
	v = v.MAdd(dt_2, a1.Add(a2))
	return p, v
}

// Symplectic Euler integrator, used for debugging.
// https://en.wikipedia.org/wiki/Semi-implicit_Euler_method
func SymEuler(p, v Vec, acc Accel, t, dt float64) (Vec, Vec) {
	v = v.MAdd(dt, acc(p, t))
	p = p.MAdd(dt, v)
	return p, v
}
