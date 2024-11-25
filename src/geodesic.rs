#[derive(Clone)]
pub struct Geodesic {
    pub t: f64,
    pub r: f64,
    pub theta: f64,
    pub phi: f64,
    pub u0: f64,
    pub u1: f64,
    pub u2: f64,
    pub u3: f64,
    pub h_step: f64,
    pub dx: f64,
    pub a: f64,
}
