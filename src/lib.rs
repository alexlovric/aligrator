pub mod adaptive_stepping;
pub mod butcher_tables;
pub mod euler;
pub mod runge_kutta;
pub mod utilities;

/// Represents a function defining an initial value problem (IVP) in the context of differential equations.
/// i.e. Second order ODE function -> xddot = f(t, x, xdot)
///
/// # Type Parameters
/// * `D`: The dimensionality of the system, indicating the size of the state and derivative vectors.
pub trait IvpFunction<const D: usize> {
    /// Computes the derivative of the state vector `x` at a given time `t`.
    ///
    /// # Arguments
    /// * `t`: The current time.
    /// * `x`: A reference to the current state vector of the system.
    /// * `xdot`: A reference to the current derivative vector of the system.
    ///
    /// # Returns
    /// * `[f64; D]`: The computed derivative vector at the given time `t`.
    fn compute(&mut self, t: &f64, x: &[f64; D], xdot: &[f64; D]) -> [f64; D];
}

/// Represents an integrator for solving initial value problems (IVPs) in differential equations.
///
/// # Type Parameters
/// * `F`: The function type that implements the `IvpFunction` trait, defining the system of differential equations.
/// * `D`: The dimensionality of the system, indicating the size of the state and derivative vectors.
pub trait Integrator<F: IvpFunction<D> + ?Sized, const D: usize> {
    /// Performs a single integration step, updating the state vector `x` and its derivative `xdot` over time `t`.
    ///
    /// # Arguments
    /// * `f`: A reference to the function defining the system of differential equations.
    /// * `t`: A mutable reference to the current time, which will be updated by the integrator.
    /// * `x`: A mutable reference to the current state vector of the system.
    /// * `xdot`: A mutable reference to the current derivative vector of the system.
    fn step(&mut self, f: &mut F, t: &mut f64, x: &mut [f64; D], xdot: &mut [f64; D]);

    /// Returns the current time step `dt` used in the integration.
    ///
    /// # Returns
    /// * `f64`: The current time step `dt`.
    fn get_dt(&self) -> f64;

    /// Sets a new time step `dt` for the integration process.
    ///
    /// # Arguments
    /// * `value`: The new time step `dt` to be set.
    fn set_dt(&mut self, value: f64);

    /// Enables or disables adaptive time-stepping in the integration process.
    ///
    /// # Arguments
    /// * `value`: A boolean value indicating whether adaptive time-stepping should be enabled (`true`) or disabled (`false`).
    fn set_adaptive(&mut self, value: bool);

    /// Returns the number of function evaluations performed per integration step.
    ///
    /// # Returns
    /// * `u16`: The number of evaluations per integration step.
    fn evaluations_per_step(&self) -> u16;
}

/// Performs the integration of a system of differential equations over a specified time interval.
///
/// # Type Parameters
/// * `I`: The integrator type, which must implement the `Integrator` trait.
/// * `F`: The function type that implements the `IvpFunction` trait, defining the system of differential equations.
/// * `D`: The dimensionality of the system, indicating the size of the state and derivative vectors.
///
/// # Arguments
/// * `integrator`: A mutable reference to the integrator that performs the numerical integration.
/// * `f`: A reference to the function that defines the system of differential equations.
/// * `x0`: The initial state vector of the system.
/// * `xdot0`: The initial derivative vector of the system.
/// * `t0`: The initial time.
/// * `tf`: The final time for the integration.
///
/// # Returns
/// * `(Vec<f64>, Vec<[f64; D]>, Vec<[f64; D]>)`: A tuple containing three vectors:
///   - A vector of time points.
///   - A vector of state vectors at each time point.
///   - A vector of derivative vectors at each time point.
pub fn integrate<I, F, const D: usize>(
    integrator: &mut I,
    f: &mut F,
    x0: [f64; D],
    xdot0: [f64; D],
    t0: f64,
    tf: f64,
) -> (Vec<f64>, Vec<[f64; D]>, Vec<[f64; D]>)
where
    I: Integrator<F, D>,
    F: IvpFunction<D>,
{
    // Initialize time and state
    let mut t = t0;
    let mut x = x0;
    let mut xdot = xdot0;

    // Vectors to store time, state values, and state derivatives
    let mut t_values = vec![t];
    let mut x_values = vec![x];
    let mut xdot_values = vec![xdot];

    while t < tf {
        // Perform the integrator step
        integrator.step(f, &mut t, &mut x, &mut xdot);

        // Store the updated time, state, and derivative
        t_values.push(t);
        x_values.push(x);
        xdot_values.push(xdot);
    }

    (t_values, x_values, xdot_values)
}
