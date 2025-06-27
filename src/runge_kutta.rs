use crate::{
    adaptive_stepping::AdaptiveDt,
    butcher_tables::{rk45_fehlberg_alfa_one_third_table, rk4_classical_table, rk89_verner_table},
    Integrator, IvpFunction,
};

#[derive(Debug)]
pub struct Rk4 {
    pub dt: f64,
    a: [f64; 4],
    b: [f64; 6],
    c: [f64; 4],
}

impl Rk4 {
    pub fn new(dt: f64) -> Self {
        let coeffs = rk4_classical_table();

        Rk4 {
            dt,
            a: coeffs.0,
            b: coeffs.1,
            c: coeffs.2,
        }
    }
}

impl<F, const D: usize> Integrator<F, D> for Rk4
where
    F: IvpFunction<D>,
{
    fn step(&mut self, f: &mut F, t: &mut f64, x: &mut [f64; D], xdot: &mut [f64; D]) {
        let dt = self.dt;
        let mut kx = [[0.; D]; 4];
        let mut k = [[0.; D]; 4];

        evaluate_rk_stages(f, t, x, xdot, &mut kx, &mut k, &dt, &self.a, &self.b);

        evaluate_rk_candidate(x, xdot, &kx, &k, &dt, &self.c);

        *t += dt;
    }

    fn get_dt(&self) -> f64 {
        self.dt
    }

    fn set_dt(&mut self, value: f64) {
        self.dt = value;
    }

    fn set_adaptive(&mut self, _: bool) {}

    fn evaluations_per_step(&self) -> u16 {
        4
    }
}

/// Rk45 Fehlberg
#[derive(Debug)]
pub struct Rk45 {
    pub dt: f64,
    adaptive_scheme: Option<AdaptiveDt>,
    adapt: bool,
    a: [f64; 6],
    b: [f64; 15],
    c: [f64; 6],
    e: [f64; 6],
}

impl Rk45 {
    pub fn new(dt: f64, mut adaptive_scheme: Option<AdaptiveDt>) -> Self {
        let coeffs = rk45_fehlberg_alfa_one_third_table();

        let mut adapt = false;
        if let Some(ref mut adaptive_scheme) = adaptive_scheme {
            adaptive_scheme.set_power(5);
            adapt = true;
        }

        Rk45 {
            dt,
            adaptive_scheme,
            adapt,
            a: coeffs.0,
            b: coeffs.1,
            c: coeffs.2,
            e: coeffs.3,
        }
    }
}

impl<F, const D: usize> Integrator<F, D> for Rk45
where
    F: IvpFunction<D>,
{
    fn step(&mut self, f: &mut F, t: &mut f64, x: &mut [f64; D], xdot: &mut [f64; D]) {
        let mut dt = self.dt;
        let mut kx = [[0.; D]; 6];
        let mut k = [[0.; D]; 6];

        match (self.adapt, &mut self.adaptive_scheme) {
            (true, Some(adaptive_scheme)) => {
                adaptive_scheme.adaptive_steps = 0;
                'adaptive_loop: loop {
                    adaptive_scheme.adaptive_steps += 1;

                    evaluate_rk_stages(f, t, x, xdot, &mut kx, &mut k, &dt, &self.a, &self.b);

                    let error = compute_rk_truncation_error(&k, &self.e);

                    if !adaptive_scheme.adapt_step_size(&mut dt, error) {
                        break 'adaptive_loop;
                    }
                    self.dt = dt;
                }
            }
            // Rk5 no adaptive time stepping
            (_, _) => evaluate_rk_stages(f, t, x, xdot, &mut kx, &mut k, &dt, &self.a, &self.b),
        }

        // Update state and time
        evaluate_rk_candidate(x, xdot, &kx, &k, &self.dt, &self.c);

        *t += self.dt;
        self.dt = dt;
    }

    fn get_dt(&self) -> f64 {
        self.dt
    }

    fn set_dt(&mut self, value: f64) {
        self.dt = value;
    }

    fn set_adaptive(&mut self, value: bool) {
        self.adapt = value;
    }

    fn evaluations_per_step(&self) -> u16 {
        match (self.adapt, &self.adaptive_scheme) {
            (true, Some(adaptive)) => 6 * adaptive.adaptive_steps,
            (_, _) => 6,
        }
    }
}

/// Otherwise Rk 89 or Dormand Prince
#[derive(Debug)]
pub struct Rk89 {
    pub dt: f64,
    adaptive_scheme: Option<AdaptiveDt>,
    adapt: bool,
    a: [f64; 16],
    b: [f64; 120],
    c: [f64; 16],
    e: [f64; 16],
}

impl Rk89 {
    pub fn new(dt: f64, mut adaptive_scheme: Option<AdaptiveDt>) -> Self {
        // let coeffs = butcher_rk5_coeffs();
        let coeffs = rk89_verner_table();

        let mut adapt = false;
        if let Some(ref mut adaptive_scheme) = adaptive_scheme {
            adaptive_scheme.set_power(9);
            adapt = true;
        }

        Rk89 {
            dt,
            adaptive_scheme,
            adapt,
            a: coeffs.0,
            b: coeffs.1,
            c: coeffs.2,
            e: coeffs.3,
        }
    }
}

impl<F, const D: usize> Integrator<F, D> for Rk89
where
    F: IvpFunction<D>,
{
    fn step(&mut self, f: &mut F, t: &mut f64, x: &mut [f64; D], xdot: &mut [f64; D]) {
        let mut dt = self.dt;
        let mut kx = [[0.; D]; 16];
        let mut k = [[0.; D]; 16];

        match (self.adapt, &mut self.adaptive_scheme) {
            (true, Some(adaptive_scheme)) => {
                adaptive_scheme.adaptive_steps = 0;
                'adaptive_loop: loop {
                    adaptive_scheme.adaptive_steps += 1;
                    evaluate_rk_stages(f, t, x, xdot, &mut kx, &mut k, &dt, &self.a, &self.b);

                    let error = compute_rk_truncation_error(&k, &self.e);

                    if !adaptive_scheme.adapt_step_size(&mut dt, error) {
                        break 'adaptive_loop;
                    }
                    self.dt = dt;
                }
            }
            // Rk9 no adaptive time stepping
            (_, _) => evaluate_rk_stages(f, t, x, xdot, &mut kx, &mut k, &dt, &self.a, &self.b),
        }

        // Update state and time
        evaluate_rk_candidate(x, xdot, &kx, &k, &self.dt, &self.c);

        *t += self.dt;
        self.dt = dt;
    }

    fn get_dt(&self) -> f64 {
        self.dt
    }

    fn set_dt(&mut self, value: f64) {
        self.dt = value;
    }

    fn set_adaptive(&mut self, value: bool) {
        self.adapt = value;
    }

    fn evaluations_per_step(&self) -> u16 {
        match (self.adapt, &self.adaptive_scheme) {
            (true, Some(adaptive)) => 16 * adaptive.adaptive_steps,
            (_, _) => 16,
        }
    }
}

// N is stages
#[allow(clippy::too_many_arguments)]
pub fn evaluate_rk_stages<const N: usize, F: IvpFunction<D>, const D: usize>(
    f: &mut F,
    t: &f64,
    x: &[f64; D],
    xdot: &[f64; D],
    kx: &mut [[f64; D]; N],
    k: &mut [[f64; D]; N],
    dt: &f64,
    a: &[f64; N],
    b: &[f64],
) {
    let mut x_cur = [0.0; D];
    let mut xdot_cur = [0.0; D];
    let mut factor: f64;
    let mut b_ind = 0;

    for st in 0..N {
        for dm in 0..D {
            kx[st][dm] = 0.;
            k[st][dm] = 0.;
            x_cur[dm] = x[dm];
            xdot_cur[dm] = xdot[dm];
        }

        for ac in 0..st {
            factor = b[b_ind] * dt;
            for dm in 0..D {
                x_cur[dm] += factor * kx[ac][dm]; // this needs to be prev xdot_cur or kx
                xdot_cur[dm] += factor * k[ac][dm];
            }
            b_ind += 1;
        }

        kx[st] = xdot_cur;
        k[st] = f.compute(&(*t + a[st] * dt), &x_cur, &xdot_cur);
    }
}

// Update candidate x, xdot
pub fn evaluate_rk_candidate<const N: usize, const D: usize>(
    x: &mut [f64; D],
    xdot: &mut [f64; D],
    kx: &[[f64; D]; N],
    k: &[[f64; D]; N],
    dt: &f64,
    c: &[f64; N],
) {
    let mut factor: f64;
    for st in 0..N {
        factor = dt * c[st];
        for dm in 0..D {
            x[dm] += factor * kx[st][dm];
            xdot[dm] += factor * k[st][dm];
        }
    }
}

pub fn compute_rk_truncation_error<const N: usize, const D: usize>(k: &[[f64; D]; N], e: &[f64; N]) -> f64 {
    let mut error = 0.0;
    for dm in 0..D {
        for st in 0..N {
            error += k[st][dm] * e[st];
        }
    }
    error.abs()
}

#[cfg(test)]
mod tests {
    use crate::{
        integrate,
        utilities::{compute_l2_error, compute_order_of_accuracy},
    };

    use super::*;

    // Simple harmonic oscillator
    fn sho_acceleration(x: f64) -> f64 {
        let omega = 1.0_f64;
        -omega.powi(2) * x
    }

    fn sho_analytical(t: f64, x: f64) -> f64 {
        let omega = 1.0;
        x * (omega * t).cos()
    }

    fn sho_initial_conditions() -> (f64, f64, f64) {
        (
            1.0,  // x0
            0.0,  // t0
            10.0, // tf
        )
    }

    struct ShoAccel;

    impl IvpFunction<1> for ShoAccel {
        fn compute(&mut self, _: &f64, x: &[f64; 1], _: &[f64; 1]) -> [f64; 1] {
            [sho_acceleration(x[0])]
        }
    }

    #[test]
    fn test_rk4_simple_harmonic_oscillator() {
        let tol = 1e-5;

        // Initial conditions
        let (x0, t0, tf) = sho_initial_conditions();
        let dt = 0.1;

        let (times, positions, _) = integrate(&mut Rk4::new(dt), &mut ShoAccel, [x0], [0.], t0, tf);

        let analytical_positions = times.iter().map(|&t| sho_analytical(t, x0)).collect::<Vec<f64>>();

        // Verify that the numerical solution is close to the analytical solution
        for (num, &anal) in positions.iter().zip(analytical_positions.iter()) {
            assert!((num[0] - anal).abs() < tol, "Numerical solution deviates from analytical solution");
        }
    }

    #[test]
    fn test_rk_order_of_accuracy() {
        // Initial conditions
        let (x0, t0, tf) = sho_initial_conditions();

        // Define different time step sizes
        let dts = vec![0.5, 0.2, 0.1];

        fn order_routine<T: Integrator<ShoAccel, 1>>(mut integrator: T, dts: &[f64], x0: f64, t0: f64, tf: f64) -> f64 {
            let mut errors: Vec<f64> = vec![];
            for &h in dts {
                integrator.set_dt(h);

                let (times, positions, _) = integrate(&mut integrator, &mut ShoAccel, [x0], [0.], t0, tf);

                let numerical_positions = positions.iter().map(|pos| pos[0]).collect::<Vec<f64>>();
                let analytical_positions = times.iter().map(|&t| sho_analytical(t, x0)).collect::<Vec<f64>>();

                // Calculate L2 error
                let error = compute_l2_error(&numerical_positions, &analytical_positions);
                errors.push(error);
            }

            // Calculate the order of accuracy
            compute_order_of_accuracy(&dts, &errors).expect("Failed to compute the order of accuracy")
        }

        // Rk4 Ensure the order of accuracy is at least 4
        let rk4_order = order_routine(Rk4::new(0.0), &dts, x0, t0, tf);
        assert!(rk4_order > 4.0, "Order of accuracy is less than 4.0: {rk4_order}");

        let rk45_order = order_routine(Rk45::new(0.0, None), &dts, x0, t0, tf);
        assert!(rk45_order > 5.0, "Order of accuracy is less than 5.0: {rk45_order}");

        let rk89_order = order_routine(Rk89::new(0.0, None), &dts, x0, t0, tf);
        assert!(rk89_order > 8.5, "Order of accuracy is less than 8.0: {rk89_order}");
    }
}
