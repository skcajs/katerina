use rand::{rngs::ThreadRng, thread_rng, Rng};

pub struct Sampler {
    rng: ThreadRng,
}

impl Sampler {
    pub fn new() -> Self {
        Sampler { rng: thread_rng() }
    }

    pub fn next(&mut self) -> f64 {
        self.rng.gen::<f64>()
    }

    pub fn next_2d(&mut self) -> (f64, f64) {
        (self.rng.gen::<f64>(), self.rng.gen::<f64>())
    }
}

#[cfg(test)]
mod tests {
    use core::f64;

    use crate::tup::Tup;

    use super::*;

    #[test]
    fn a_sampler() {
        let mut s = Sampler::new();
        let n = s.next();
        assert!(n >= 0.0 && n <= 1.0);
    }

    #[test]
    fn a_2d_sampler() {
        let mut s = Sampler::new();
        let (n1, n2) = s.next_2d();
        assert!(n1 >= 0.0 && n1 <= 1.0);
        assert!(n2 >= 0.0 && n2 <= 1.0);
    }

    #[test]
    fn scalar_sample() {
        let mut s = Sampler::new();
        let mut count = 0.;
        for _ in 0..1000000 {
            let theta = s.next() * f64::consts::PI;
            let phi = s.next() * 2.0 * f64::consts::PI;

            let tup: Tup = Tup(1., theta, phi);
            let tup_cart = tup.spherical_to_cartesian();
            count += tup_cart.0;
        }

        println!("Count: {}", count);
    }
}
