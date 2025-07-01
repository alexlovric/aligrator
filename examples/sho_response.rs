#![allow(unused_imports)]

use std::fs::File;
use std::io::{BufWriter, Result, Write};

use aligrator::adaptive_stepping::AdaptiveDt;
use aligrator::runge_kutta::{Rk4, Rk45, Rk89};
use aligrator::{integrate, Integrator, IvpFunction};

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

fn main() -> Result<()> {
    // Initial conditions
    let x0 = [1.0];
    let xdot0 = [0.0];
    let t0 = 0.0;
    let tf = 15.0;
    let dt = 0.1;

    // Run the Forward Euler method
    // let mut integrator = ForwardEuler::new(dt);
    // let mut integrator = Rk4::new(dt);
    // let mut integrator = Rk45::new(dt, Some(AdaptiveDt::new(Some(1e-4), None, None)));
    let mut integrator = Rk89::new(dt, None);
    // let mut integrator = Rk89::new(dt, Some(AdaptiveDt::new(Some(1e-6), None, None)));

    let (times, positions, _) = integrate(
        &mut integrator,
        &mut SimpleHarmonicOscillator { omega: 1.0 },
        x0,
        xdot0,
        t0,
        tf,
    );

    // Write out data
    let file = File::create("sho_response.csv")?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "time,position.x,position.y,position.z")?;

    for (i, &t) in times.iter().enumerate() {
        writeln!(writer, "{},{}", t, positions[i][0])?;
    }

    writer.flush()?;

    println!("Response written to 'sho_response.csv'");

    Ok(())
}
