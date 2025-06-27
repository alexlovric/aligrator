use crate::{Integrator, IvpFunction};

/// Perform the Forward Euler method to solve the ODE given by a function.
#[derive(Debug)]
pub struct ForwardEuler {
    dt: f64,
}

impl ForwardEuler {
    pub fn new(dt: f64) -> Self {
        ForwardEuler { dt }
    }
}

impl<F, const D: usize> Integrator<F, D> for ForwardEuler
where
    F: IvpFunction<D>,
{
    fn step(&mut self, f: &mut F, t: &mut f64, x: &mut [f64; D], xdot: &mut [f64; D]) {
        let a = f.compute(t, x, xdot); // xddot_n = f(t_n, x_n, xdot_n)

        // Update position and velocity using the factor
        for i in 0..D {
            x[i] += xdot[i] * self.dt; // x_(n+1) = x_n + xdot_n * dt
            xdot[i] += a[i] * self.dt; // xdot_(n+1) = xdot_n + xddot_n * dt
        }

        // Can adjust dt if necessary
        *t += self.dt; // t_(n+1) = t_n + dt
    }

    fn get_dt(&self) -> f64 {
        self.dt
    }

    fn set_dt(&mut self, value: f64) {
        self.dt = value;
    }

    fn set_adaptive(&mut self, _: bool) {}

    fn evaluations_per_step(&self) -> u16 {
        1
    }
}
