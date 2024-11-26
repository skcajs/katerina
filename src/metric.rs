use crate::{geodesic::Geodesic, tup::Tup};

#[derive(Clone)]
pub struct Metric {
    pub s: Tup,
    pub a: f64,
    pub rs: f64,
}

impl Metric {
    pub fn new(rs: f64, s: Tup, a: f64) -> Self {
        Metric { rs: rs * 0.6, s, a } // applying a scale factor to the Schwarzschild radius
    }

    #[allow(dead_code)]
    pub fn minkowski(&self, geo: &Geodesic, h: f64) -> Geodesic {
        let current_point = geo.ray.o + geo.ray.d * h;

        Geodesic::ray(current_point, geo.ray.d, self.clone())
    }

    #[allow(dead_code)]
    pub fn rkf45(&self, geo: &Geodesic, h: f64, tol: f64) -> (Geodesic, f64) {
        let mut h = h;
        let previous_x = geo.ray.o - self.s;
        let previous_p = geo.ray.d;

        loop {
            let (k1x, k1p) = self.schwarzschild(previous_p, previous_x);
            let (k2x, k2p) =
                self.schwarzschild(previous_p + k1p * 0.25 * h, previous_x + k1x * 0.25 * h);
            let (k3x, k3p) = self.schwarzschild(
                previous_p + (k1p * (3. / 32.) + k2p * (9. / 32.)) * h,
                previous_x + (k1x * (3. / 32.) + k2x * (9. / 32.)) * h,
            );
            let (k4x, k4p) = self.schwarzschild(
                previous_p
                    + (k1p * (1932. / 2197.) - k2p * (7200. / 2197.) + k3p * (7296. / 2197.)) * h,
                previous_x
                    + (k1x * (1932. / 2197.) - k2x * (7200. / 2197.) + k3x * (7296. / 2197.)) * h,
            );
            let (k5x, k5p) = self.schwarzschild(
                previous_p
                    + (k1p * (439. / 216.) - k2p * 8. + k3p * (3680. / 513.)
                        - k4p * (845. / 4104.))
                        * h,
                previous_x
                    + (k1x * (439. / 216.) - k2x * 8. + k3x * (3680. / 513.)
                        - k4x * (845. / 4104.))
                        * h,
            );
            let (k6x, k6p) = self.schwarzschild(
                previous_p
                    - (k1p * (8. / 27.) + k2p * 2. - k3p * (3544. / 2565.) + k4p * (1859. / 4104.)
                        - k5p * (11. / 40.))
                        * h,
                previous_x
                    - (k1x * (8. / 27.) + k2x * 2. - k3x * (3544. / 2565.) + k4x * (1859. / 4104.)
                        - k5x * (11. / 40.))
                        * h,
            );

            // 4th-order solution
            let x4 = previous_x
                + (k1x * (25. / 216.) * h + k3x * (1408. / 2565.) * h + k4x * (2197. / 4104.) * h
                    - k5x * (1. / 5.) * h);
            let p4 = previous_p
                + (k1p * (25. / 216.) * h + k3p * (1408. / 2565.) * h + k4p * (2197. / 4104.) * h
                    - k5p * (1. / 5.) * h);

            // 5th-order solution
            let x5 = previous_x
                + (k1x * (16. / 135.) + k3x * (6656. / 12825.) + k4x * (28561. / 56430.)
                    - k5x * (9. / 50.)
                    + k6x * (2. / 55.))
                    * h;
            let p5 = previous_p
                + (k1p * (16. / 135.) + k3p * (6656. / 12825.) + k4p * (28561. / 56430.)
                    - k5p * (9. / 50.)
                    + k6p * (2. / 55.))
                    * h;

            // Error estimation
            let error_x = (x4 - x5).len();
            let error_p = (p4 - p5).len();
            let error = error_x.max(error_p);

            // Adjust step size based on the error
            if error < tol {
                return (Geodesic::ray(x5 + self.s, p5.norm(), self.clone()), h);
            }

            // Reduce step size for better accuracy
            h *= 0.9 * (tol / error).powf(0.2);
        }
    }

    #[allow(dead_code)]
    pub fn rk4(&self, geo: &Geodesic, h: f64) -> Geodesic {
        // println!("h: {}", h);
        let previous_x = geo.ray.o - self.s;
        let previous_p = geo.ray.d;

        let (k1x, k1p) = self.schwarzschild(previous_p, previous_x);
        let (k2x, k2p) = self.schwarzschild(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);
        let (k3x, k3p) = self.schwarzschild(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);
        let (k4x, k4p) = self.schwarzschild(previous_p + k3p * h, previous_x + k3x * h);

        let current_point = previous_x + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
        let current_momentum = previous_p + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

        Geodesic::ray(
            current_point + self.s,
            current_momentum.norm(),
            self.clone(),
        )
    }

    #[allow(dead_code)]
    pub fn euler(&self, geo: &Geodesic, h: f64) -> Geodesic {
        let pos = geo.ray.o;
        let dir = geo.ray.d.norm();

        let r: f64 = (pos - self.s).len();

        if r < 1. {
            return Geodesic::ray(pos, dir, self.clone());
        }

        let h2 = pos.cross(dir).len().powi(2);

        let new_pos = pos + dir * h;
        let new_dir = dir + (pos * -1.5 * h2 * (1. / r.powi(2).powf(2.5) * h));

        Geodesic::ray(new_pos, new_dir, self.clone())
    }

    fn schwarzschild(&self, p: Tup, x: Tup) -> (Tup, Tup) {
        let r: f64 = x.len();
        if r < self.rs / 4. {
            return (p, x);
        }
        let r_adjusted = r * (1. + (self.rs / (4. * r))).powi(2);
        let a: f64 = 1. + (self.rs / (4. * r));
        let b: f64 = 1. - (self.rs / (4. * r));
        let fact_x: f64 = (b * b) / a.powi(6);
        let fact_p1: f64 = (b * b) / a.powi(7);
        let fact_p2: f64 = 1.0 / (b * a);
        (
            p * fact_x,
            x * (-1. / (2. * r_adjusted.powi(3)))
                * (((p.0.powi(2) + p.1.powi(2) + p.2.powi(2)) * fact_p1) + fact_p2)
                * self.rs,
        )
    }

    pub fn transform_point(&self, light_pos: Tup, observer_pos: Tup) -> Tup {
        // Shift positions relative to the black hole's position
        let light_rel = light_pos - self.s;
        let observer_rel = observer_pos - self.s;

        // Compute radial distance in the black hole's local coordinate system
        let r = (light_rel - observer_rel).len();
        let factor = 1. + (self.rs / (4. * r)).powi(4);

        // Adjust radial distance
        let adjusted_r = r * factor;

        // Transform point in the black hole's local coordinates
        let transformed_rel = observer_rel + (light_rel - observer_rel).norm() * adjusted_r;

        // Convert back to world coordinates
        transformed_rel + self.s
    }
}
