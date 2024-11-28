use crate::tup::Tup;

#[derive(Clone, Debug)]
pub struct Metric {
    pub s: Tup,
    pub a: f64,
    pub rs: f64,
    pub initial_step: f64,
}

impl Metric {
    pub fn new(s: Tup, a: f64, m: f64) -> Self {
        Metric {
            rs: 2. * m * 0.6,
            s,
            a,
            initial_step: 0.,
        } // applying a scale factor to the Schwarzschild radius
    }
}
