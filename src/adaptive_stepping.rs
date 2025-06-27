#[derive(Debug)]
pub struct AdaptiveDt {
    pub tol: f64,
    pub min: f64,
    pub max: Option<f64>,
    power_increase: f64,
    power_decrease: f64,
    safety: f64,
    pub adaptive_steps: u16,
}

impl AdaptiveDt {
    pub fn new(tol: Option<f64>, min: Option<f64>, max: Option<f64>) -> Self {
        AdaptiveDt {
            tol: tol.unwrap_or(1e-5),
            min: min.unwrap_or(1e-20), // small enough?
            max,
            power_increase: 1.0,
            power_decrease: 1.0,
            safety: 0.90,
            adaptive_steps: 0,
        }
    }

    pub fn set_power(&mut self, order: u8) {
        self.power_increase = 1.0 / (order as f64);
        self.power_decrease = 1.0 / ((order - 1) as f64);
    }

    pub fn adapt_step_size(&mut self, dt: &mut f64, error: f64) -> bool {
        let ratio = if error < 1e-22 {
            println!("🟠 Truncation error is {error:.4e}, hence dividing by zero! Setting as 1e-20");
            self.tol / 1e-20 // some default
        } else {
            self.tol / error
        };

        let mut flag = false;

        let factor = if error > self.tol {
            if self.adaptive_steps > 9 {
                println!("🟠 Adaptive iteration {}, exceeds limit!", self.adaptive_steps);
            } else {
                flag = true;
            }

            ratio.powf(self.power_decrease)
        } else {
            ratio.powf(self.power_increase)
        };

        *dt *= self.safety * factor;
        *dt = if let Some(max) = self.max {
            dt.clamp(self.min, max)
        } else {
            dt.max(self.min)
        };

        flag
    }
}

impl Default for AdaptiveDt {
    fn default() -> Self {
        AdaptiveDt::new(Some(1e-9), None, None)
    }
}
