#![allow(unused_imports)]

use std::fs::File;
use std::io::{self, BufWriter, Write};

use aligrator::runge_kutta::Rk89;
use aligrator::utilities::compute_l2_error;
use aligrator::{Integrator, IvpFunction, integrate};

/// Simple Harmonic Oscillator (SHO) implementation.
///
/// This implements the equation of motion for a simple harmonic oscillator:
///     x''(t) + ω²x(t) = 0
///
/// Where:
/// - x(t) is the displacement from equilibrium at time t
/// - ω is the angular frequency of oscillation (ω = 2πf, where f is the frequency)
/// - x''(t) is the second time derivative of x (acceleration)
struct SimpleHarmonicOscillator {
    omega: f64,
}

impl IvpFunction<1> for SimpleHarmonicOscillator {
    fn compute(&mut self, _: &f64, x: &[f64; 1], _: &[f64; 1]) -> [f64; 1] {
        [-self.omega.powi(2) * x[0]]
    }
}

fn main() -> io::Result<()> {
    // SHO parameters
    let omega = 1.0; // angular frequency

    // Initial conditions
    let x0 = [1.0];
    let xdot0 = [0.0];
    let t0 = 0.0;
    let tf = 10.0;

    // Different delta t values to test
    let delta_ts = vec![0.5, 0.2, 0.1];

    // let mut integrator = ForwardEuler::new(0.0);
    // let mut integrator = Rk4::new(0.0);
    // let mut integrator = Rk45::new(0.0, Some(AdaptiveDt::new(Some(1e-5), Some(0.1), Some(0.5))));
    // let mut integrator = Rk45::new(0.0, None);
    let mut integrator = Rk89::new(0.0, None);

    let mut results = vec![];

    for &h in &delta_ts {
        integrator.dt = h;
        let (times, positions, _) = integrate(&mut integrator, &mut SimpleHarmonicOscillator { omega: 1.0 }, x0, xdot0, t0, tf);

        let numerical_positions: Vec<f64> = positions.iter().map(|pos| pos[0]).collect();

        // Calculate analytical positions
        let analytical_positions: Vec<f64> =
            times.iter().map(|&t| x0[0] * (omega * t).cos()).collect();

        // Calculate L2 error
        let error = compute_l2_error(&numerical_positions, &analytical_positions);
        results.push((h, error));
    }

    let file = File::create("sho_errors.csv")?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "delta_t,l2_error")?;
    for (delta_t, error) in results {
        writeln!(writer, "{},{}", delta_t, error)?;
    }
    writer.flush()?;

    println!("Errors written to 'sho_errors.csv'");

    Ok(())
}
