use nalgebra::SVector;

use crate::{geodesic::Geodesic, ray::Ray, tup::Tup};

#[derive(Clone, Debug)]
pub struct Metric {
    pub s: Tup,
    pub a: f64,
    pub rs: f64,
}

impl Metric {
    pub fn new(s: Tup, a: f64, m: f64) -> Self {
        Metric {
            rs: 2. * m * 0.6,
            s,
            a,
        } // applying a scale factor to the Schwarzschild radius
    }

    #[allow(dead_code)]
    pub fn minkowski(&self, geo: &Geodesic, h: f64) -> Geodesic {
        let current_point = geo.ray.o + geo.ray.d * h;

        Geodesic::ray(current_point, geo.ray.d, self.clone())
    }

    #[allow(dead_code)]
    pub fn rk4_kerr(&self, geo: &Geodesic, h: f64) -> Geodesic {
        let incoming: SVector<f64, 8> = SVector::<f64, 8>::from_vec(vec![
            geo.t, geo.r, geo.theta, geo.phi, geo.u0, geo.u1, geo.u2, geo.u3,
        ]);

        let k1: SVector<f64, 8> = self.kerr(incoming);
        let k2: SVector<f64, 8> = self.kerr(incoming + k1 * 0.5 * h);
        let k3: SVector<f64, 8> = self.kerr(incoming + k2 * 0.5 * h);
        let k4: SVector<f64, 8> = self.kerr(incoming + k3 * h);

        let current_point: SVector<f64, 8> = incoming + (k1 + k2 * 2. + k3 * 2. + k4) * (h / 6.);

        let point = Tup(current_point[1], current_point[2], current_point[3])
            .boyer_lindquist_to_cartesian(geo.m.a)
            + self.s;

        let dir = (point - geo.ray.o).norm();

        Geodesic::ray(point, dir, self.clone())
    }

    #[allow(dead_code)]
    pub fn rkf45_kerr(&self, geo: &Geodesic, h: &mut f64) -> Geodesic {
        let tol = 1e-6;
        let incoming: SVector<f64, 8> = SVector::<f64, 8>::from_vec(vec![
            geo.t, geo.r, geo.theta, geo.phi, geo.u0, geo.u1, geo.u2, geo.u3,
        ]);

        loop {
            let k1: SVector<f64, 8> = self.kerr(incoming);
            let k2: SVector<f64, 8> = self.kerr(incoming + k1 * 0.25 * *h);
            let k3: SVector<f64, 8> =
                self.kerr(incoming + (k1 * (3. / 32.) + k2 * (9. / 32.)) * *h);
            let k4: SVector<f64, 8> = self.kerr(
                incoming
                    + (k1 * (1932. / 2197.) - k2 * (7200. / 2197.) + k3 * (7296. / 2197.)) * *h,
            );
            let k5: SVector<f64, 8> = self.kerr(
                incoming
                    + (k1 * (439. / 216.) - k2 * 8. + k3 * (3680. / 513.) - k4 * (845. / 4104.))
                        * *h,
            );
            let k6: SVector<f64, 8> = self.kerr(
                incoming
                    - (k1 * (8. / 27.) + k2 * 2. - k3 * (3544. / 2565.) + k4 * (1859. / 4104.)
                        - k5 * (11. / 40.))
                        * *h,
            );

            // 4th-order solution
            let x4: SVector<f64, 8> = incoming
                + (k1 * (25. / 216.) * *h + k3 * (1408. / 2565.) * *h + k4 * (2197. / 4104.) * *h
                    - k5 * (1. / 5.) * *h);

            // 5th-order solution
            let x5: SVector<f64, 8> = incoming
                + (k1 * (16. / 135.) + k3 * (6656. / 12825.) + k4 * (28561. / 56430.)
                    - k5 * (9. / 50.)
                    + k6 * (2. / 55.))
                    * *h;

            // Error estimation
            let error_val: SVector<f64, 8> = x4 - x5;
            let error = error_val.norm();

            // Adjust step size based on the error
            if error < tol {
                let point = Tup(x5[1], x5[2], x5[3]).boyer_lindquist_to_cartesian(geo.m.a) + self.s;
                let dir = (point - geo.ray.o).norm();
                return Geodesic::ray(point, dir, self.clone());
            }

            // Reduce step size for better accuracy
            *h *= 0.9 * (tol / error).powf(0.2);
        }
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

    #[allow(dead_code)]
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

    #[allow(dead_code)]
    fn kerr(&self, incoming: SVector<f64, 8>) -> SVector<f64, 8> {
        let a = self.a;
        // let t = incoming[0];
        let r = incoming[1];
        let theta = incoming[2];
        // let phi = incoming[3];
        let u0 = incoming[4];
        let u1 = incoming[5];
        let u2 = incoming[6];
        let u3 = incoming[7];

        let r2 = r * r;
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a2 * a2;

        let cos_theta = f64::cos(theta);
        let sin_theta = f64::sin(theta);
        let cos_2_theta = cos_theta * cos_theta;
        let sin_2_theta = sin_theta * sin_theta;
        let sin_3_theta = sin_2_theta * sin_theta;
        let s = a2 * cos_2_theta + r2;
        let s2 = s * s;

        let u0_dot = -2.0
            * u0
            * u1
            * (-a * r * (4.0 * a * r2 * sin_2_theta / s2 - 2.0 * a * sin_2_theta / s)
                / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                    + a4 * cos_2_theta
                    + a2 * r2 * cos_2_theta
                    + a2 * r2
                    + r4)
                + 0.5
                    * (-4.0 * r2 / s2 + 2.0 / s)
                    * (-2.0 * a2 * r * sin_2_theta
                        - a4 * cos_2_theta
                        - a2 * r2 * cos_2_theta
                        - a2 * r2
                        - r4)
                    / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_2_theta
                        + a2 * r2 * cos_2_theta
                        + a2 * r2
                        + r4))
            - 2.0
                * u0
                * u2
                * (2.0
                    * a2
                    * r
                    * (-2.0 * a2 * r * sin_2_theta
                        - a4 * cos_2_theta
                        - a2 * r2 * cos_2_theta
                        - a2 * r2
                        - r4)
                    * sin_theta
                    * cos_theta
                    / (s2
                        * (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_2_theta
                            + a2 * r2 * cos_2_theta
                            + a2 * r2
                            + r4))
                    - a * r
                        * (-4.0 * a3 * r * sin_3_theta * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / s)
                        / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_2_theta
                            + a2 * r2 * cos_2_theta
                            + a2 * r2
                            + r4))
            - 2.0
                * u1
                * u3
                * (-a
                    * r
                    * (-4.0 * a2 * r2 * sin_2_theta / s2 + 2.0 * a2 * sin_2_theta / s + 2.0 * r)
                    * sin_2_theta
                    / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_2_theta
                        + a2 * r2 * cos_2_theta
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (4.0 * a * r2 * sin_2_theta / s2 - 2.0 * a * sin_2_theta / s)
                        * (-2.0 * a2 * r * sin_2_theta
                            - a4 * cos_2_theta
                            - a2 * r2 * cos_2_theta
                            - a2 * r2
                            - r4)
                        / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_2_theta
                            + a2 * r2 * cos_2_theta
                            + a2 * r2
                            + r4))
            - 2.0
                * u2
                * u3
                * (-a
                    * r
                    * ((4.0 * a4 * r * sin_3_theta * cos_theta / s2
                        + 4.0 * a2 * r * sin_theta * cos_theta / s)
                        * sin_2_theta
                        + 2.0
                            * (2.0 * a2 * r * sin_2_theta / s + a2 + r2)
                            * sin_theta
                            * cos_theta)
                    / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_2_theta
                        + a2 * r2 * cos_2_theta
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (-4.0 * a3 * r * sin_3_theta * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / s)
                        * (-2.0 * a2 * r * sin_2_theta
                            - a4 * cos_2_theta
                            - a2 * r2 * cos_2_theta
                            - a2 * r2
                            - r4)
                        / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_2_theta
                            + a2 * r2 * cos_2_theta
                            + a2 * r2
                            + r4));
        let u1_dot = 2.0 * a2 * u1 * u2 * sin_theta * cos_theta / s
            + r * u2 * u2 * (-2.0 * r + a2 + r2) / s
            - 0.5 * u0 * u0 * (4.0 * r2 / s2 - 2.0 / s) * (-2.0 * r + a2 + r2) / s
            - u0 * u3
                * (-4.0 * a * r2 * sin_2_theta / s2 + 2.0 * a * sin_2_theta / s)
                * (-2.0 * r + a2 + r2)
                / s
            - 0.5
                * u1
                * u1
                * (2.0 * r / (-2.0 * r + a2 + r2)
                    + (2.0 - 2.0 * r) * s / ((-2.0 * r + a2 + r2) * (-2.0 * r + a2 + r2)))
                * (-2.0 * r + a2 + r2)
                / s
            + 0.5
                * u3
                * u3
                * (-2.0 * r + a2 + r2)
                * (-4.0 * a2 * r2 * sin_2_theta / s2 + 2.0 * a2 * sin_2_theta / s + 2.0 * r)
                * sin_2_theta
                / s;
        let u2_dot = 2.0 * a2 * r * u0 * u0 * sin_theta * cos_theta / (s2 * s)
            - a2 * u1 * u1 * sin_theta * cos_theta / ((s) * (-2.0 * r + a2 + r2))
            + a2 * u2 * u2 * sin_theta * cos_theta / (s)
            - 2.0 * r * u1 * u2 / (s)
            - u0 * u3
                * (4.0 * a3 * r * sin_3_theta * cos_theta / s2
                    + 4.0 * a * r * sin_theta * cos_theta / (s))
                / (s)
            - 0.5
                * u3
                * u3
                * (-(4.0 * a4 * r * sin_3_theta * cos_theta / s2
                    + 4.0 * a2 * r * sin_theta * cos_theta / (s))
                    * sin_2_theta
                    - 2.0 * (2.0 * a2 * r * sin_2_theta / (s) + a2 + r2) * sin_theta * cos_theta)
                / (s);
        let u3_dot = -2.0
            * u0
            * u1
            * (-a * r * (-4.0 * r2 / s2 + 2.0 / (s))
                / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                    + a4 * cos_2_theta
                    + a2 * r2 * cos_2_theta
                    + a2 * r2
                    + r4)
                + 0.5
                    * (4.0 * a * r2 * sin_2_theta / s2 - 2.0 * a * sin_2_theta / (s))
                    * (-2.0 * r + s)
                    / (2.0 * a2 * r * sin_2_theta * sin_2_theta
                        - 2.0 * a2 * r * sin_2_theta
                        - 2.0 * r3 * sin_2_theta
                        + a4 * sin_2_theta * cos_2_theta
                        + a2 * r2 * sin_2_theta * cos_2_theta
                        + a2 * r2 * sin_2_theta
                        + r4 * sin_2_theta))
            - 2.0
                * u0
                * u2
                * (-4.0 * a3 * r2 * sin_theta * cos_theta
                    / (s2
                        * (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_2_theta
                            + a2 * r2 * cos_2_theta
                            + a2 * r2
                            + r4))
                    + 0.5
                        * (-4.0 * a3 * r * sin_3_theta * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / (s))
                        * (-2.0 * r + s)
                        / (2.0 * a2 * r * sin_2_theta * sin_2_theta
                            - 2.0 * a2 * r * sin_2_theta
                            - 2.0 * r3 * sin_2_theta
                            + a4 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta
                            + r4 * sin_2_theta))
            - 2.0
                * u1
                * u3
                * (-a * r * (4.0 * a * r2 * sin_2_theta / s2 - 2.0 * a * sin_2_theta / (s))
                    / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_2_theta
                        + a2 * r2 * cos_2_theta
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (-2.0 * r + s)
                        * (-4.0 * a2 * r2 * sin_2_theta / s2
                            + 2.0 * a2 * sin_2_theta / (s)
                            + 2.0 * r)
                        * sin_2_theta
                        / (2.0 * a2 * r * sin_2_theta * sin_2_theta
                            - 2.0 * a2 * r * sin_2_theta
                            - 2.0 * r3 * sin_2_theta
                            + a4 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta
                            + r4 * sin_2_theta))
            - 2.0
                * u2
                * u3
                * (-a
                    * r
                    * (-4.0 * a3 * r * sin_3_theta * cos_theta / s2
                        - 4.0 * a * r * sin_theta * cos_theta / (s))
                    / (2.0 * a2 * r * sin_2_theta - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_2_theta
                        + a2 * r2 * cos_2_theta
                        + a2 * r2
                        + r4)
                    + 0.5
                        * ((4.0 * a4 * r * sin_3_theta * cos_theta / s2
                            + 4.0 * a2 * r * sin_theta * cos_theta / (s))
                            * sin_2_theta
                            + 2.0
                                * (2.0 * a2 * r * sin_2_theta / (s) + a2 + r2)
                                * sin_theta
                                * cos_theta)
                        * (-2.0 * r + s)
                        / (2.0 * a2 * r * sin_2_theta * sin_2_theta
                            - 2.0 * a2 * r * sin_2_theta
                            - 2.0 * r3 * sin_2_theta
                            + a4 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta * cos_2_theta
                            + a2 * r2 * sin_2_theta
                            + r4 * sin_2_theta));

        SVector::<f64, 8>::from_vec(vec![u0, u1, u2, u3, u0_dot, u1_dot, u2_dot, u3_dot])
    }

    #[allow(dead_code)]
    fn derivatives(&self, incoming: SVector<f64, 8>) -> SVector<f64, 8> {
        let a = self.a;
        let t = incoming[0];
        let r = incoming[1];
        let theta = incoming[2];
        let phi = incoming[3];
        let u0 = incoming[4];
        let u1 = incoming[5];
        let u2 = incoming[6];
        let u3 = incoming[7];

        let r2 = r * r;
        let r3 = r2 * r;
        let r4 = r2 * r2;
        let a2 = a * a;
        let a3 = a2 * a;
        let a4 = a2 * a2;
        let cos_theta = theta.cos();
        let sin_theta = theta.sin();
        let cos_theta2 = cos_theta * cos_theta;
        let sin_theta2 = sin_theta * sin_theta;
        let sin_theta3 = sin_theta2 * sin_theta;
        let s = a2 * cos_theta2 + r2;
        let s2 = s * s;

        let u0_dot = -2.0
            * u0
            * u1
            * (-a * r * (4.0 * a * r2 * sin_theta2 / s2 - 2.0 * a * sin_theta2 / s)
                / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                    + a4 * cos_theta2
                    + a2 * r2 * cos_theta2
                    + a2 * r2
                    + r4)
                + 0.5
                    * (-4.0 * r2 / s2 + 2.0 / s)
                    * (-2.0 * a2 * r * sin_theta2
                        - a4 * cos_theta2
                        - a2 * r2 * cos_theta2
                        - a2 * r2
                        - r4)
                    / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_theta2
                        + a2 * r2 * cos_theta2
                        + a2 * r2
                        + r4))
            - 2.0
                * u0
                * u2
                * (2.0
                    * a2
                    * r
                    * (-2.0 * a2 * r * sin_theta2
                        - a4 * cos_theta2
                        - a2 * r2 * cos_theta2
                        - a2 * r2
                        - r4)
                    * sin_theta
                    * cos_theta
                    / (s2
                        * (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_theta2
                            + a2 * r2 * cos_theta2
                            + a2 * r2
                            + r4))
                    - a * r
                        * (-4.0 * a3 * r * sin_theta3 * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / s)
                        / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_theta2
                            + a2 * r2 * cos_theta2
                            + a2 * r2
                            + r4))
            - 2.0
                * u1
                * u3
                * (-a
                    * r
                    * (-4.0 * a2 * r2 * sin_theta2 / s2 + 2.0 * a2 * sin_theta2 / s + 2.0 * r)
                    * sin_theta2
                    / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_theta2
                        + a2 * r2 * cos_theta2
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (4.0 * a * r2 * sin_theta2 / s2 - 2.0 * a * sin_theta2 / s)
                        * (-2.0 * a2 * r * sin_theta2
                            - a4 * cos_theta2
                            - a2 * r2 * cos_theta2
                            - a2 * r2
                            - r4)
                        / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_theta2
                            + a2 * r2 * cos_theta2
                            + a2 * r2
                            + r4))
            - 2.0
                * u2
                * u3
                * (-a
                    * r
                    * ((4.0 * a4 * r * sin_theta3 * cos_theta / s2
                        + 4.0 * a2 * r * sin_theta * cos_theta / s)
                        * sin_theta2
                        + 2.0 * (2.0 * a2 * r * sin_theta2 / s + a2 + r2) * sin_theta * cos_theta)
                    / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_theta2
                        + a2 * r2 * cos_theta2
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (-4.0 * a3 * r * sin_theta3 * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / s)
                        * (-2.0 * a2 * r * sin_theta2
                            - a4 * cos_theta2
                            - a2 * r2 * cos_theta2
                            - a2 * r2
                            - r4)
                        / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_theta2
                            + a2 * r2 * cos_theta2
                            + a2 * r2
                            + r4));
        let u1_dot = 2.0 * a2 * u1 * u2 * sin_theta * cos_theta / s
            + r * u2 * u2 * (-2.0 * r + a2 + r2) / s
            - 0.5 * u0 * u0 * (4.0 * r2 / s2 - 2.0 / s) * (-2.0 * r + a2 + r2) / s
            - u0 * u3
                * (-4.0 * a * r2 * sin_theta2 / s2 + 2.0 * a * sin_theta2 / s)
                * (-2.0 * r + a2 + r2)
                / s
            - 0.5
                * u1
                * u1
                * (2.0 * r / (-2.0 * r + a2 + r2)
                    + (2.0 - 2.0 * r) * s / ((-2.0 * r + a2 + r2) * (-2.0 * r + a2 + r2)))
                * (-2.0 * r + a2 + r2)
                / s
            + 0.5
                * u3
                * u3
                * (-2.0 * r + a2 + r2)
                * (-4.0 * a2 * r2 * sin_theta2 / s2 + 2.0 * a2 * sin_theta2 / s + 2.0 * r)
                * sin_theta2
                / s;
        let u2_dot = 2.0 * a2 * r * u0 * u0 * sin_theta * cos_theta / (s2 * s)
            - a2 * u1 * u1 * sin_theta * cos_theta / ((s) * (-2.0 * r + a2 + r2))
            + a2 * u2 * u2 * sin_theta * cos_theta / (s)
            - 2.0 * r * u1 * u2 / (s)
            - u0 * u3
                * (4.0 * a3 * r * sin_theta3 * cos_theta / s2
                    + 4.0 * a * r * sin_theta * cos_theta / (s))
                / (s)
            - 0.5
                * u3
                * u3
                * (-(4.0 * a4 * r * sin_theta3 * cos_theta / s2
                    + 4.0 * a2 * r * sin_theta * cos_theta / (s))
                    * sin_theta2
                    - 2.0 * (2.0 * a2 * r * sin_theta2 / (s) + a2 + r2) * sin_theta * cos_theta)
                / (s);
        let u3_dot = -2.0
            * u0
            * u1
            * (-a * r * (-4.0 * r2 / s2 + 2.0 / (s))
                / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                    + a4 * cos_theta2
                    + a2 * r2 * cos_theta2
                    + a2 * r2
                    + r4)
                + 0.5
                    * (4.0 * a * r2 * sin_theta2 / s2 - 2.0 * a * sin_theta2 / (s))
                    * (-2.0 * r + s)
                    / (2.0 * a2 * r * sin_theta2 * sin_theta2
                        - 2.0 * a2 * r * sin_theta2
                        - 2.0 * r3 * sin_theta2
                        + a4 * sin_theta2 * cos_theta2
                        + a2 * r2 * sin_theta2 * cos_theta2
                        + a2 * r2 * sin_theta2
                        + r4 * sin_theta2))
            - 2.0
                * u0
                * u2
                * (-4.0 * a3 * r2 * sin_theta * cos_theta
                    / (s2
                        * (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                            + a4 * cos_theta2
                            + a2 * r2 * cos_theta2
                            + a2 * r2
                            + r4))
                    + 0.5
                        * (-4.0 * a3 * r * sin_theta3 * cos_theta / s2
                            - 4.0 * a * r * sin_theta * cos_theta / (s))
                        * (-2.0 * r + s)
                        / (2.0 * a2 * r * sin_theta2 * sin_theta2
                            - 2.0 * a2 * r * sin_theta2
                            - 2.0 * r3 * sin_theta2
                            + a4 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2
                            + r4 * sin_theta2))
            - 2.0
                * u1
                * u3
                * (-a * r * (4.0 * a * r2 * sin_theta2 / s2 - 2.0 * a * sin_theta2 / (s))
                    / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_theta2
                        + a2 * r2 * cos_theta2
                        + a2 * r2
                        + r4)
                    + 0.5
                        * (-2.0 * r + s)
                        * (-4.0 * a2 * r2 * sin_theta2 / s2
                            + 2.0 * a2 * sin_theta2 / (s)
                            + 2.0 * r)
                        * sin_theta2
                        / (2.0 * a2 * r * sin_theta2 * sin_theta2
                            - 2.0 * a2 * r * sin_theta2
                            - 2.0 * r3 * sin_theta2
                            + a4 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2
                            + r4 * sin_theta2))
            - 2.0
                * u2
                * u3
                * (-a
                    * r
                    * (-4.0 * a3 * r * sin_theta3 * cos_theta / s2
                        - 4.0 * a * r * sin_theta * cos_theta / (s))
                    / (2.0 * a2 * r * sin_theta2 - 2.0 * a2 * r - 2.0 * r3
                        + a4 * cos_theta2
                        + a2 * r2 * cos_theta2
                        + a2 * r2
                        + r4)
                    + 0.5
                        * ((4.0 * a4 * r * sin_theta3 * cos_theta / s2
                            + 4.0 * a2 * r * sin_theta * cos_theta / (s))
                            * sin_theta2
                            + 2.0
                                * (2.0 * a2 * r * sin_theta2 / (s) + a2 + r2)
                                * sin_theta
                                * cos_theta)
                        * (-2.0 * r + s)
                        / (2.0 * a2 * r * sin_theta2 * sin_theta2
                            - 2.0 * a2 * r * sin_theta2
                            - 2.0 * r3 * sin_theta2
                            + a4 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2 * cos_theta2
                            + a2 * r2 * sin_theta2
                            + r4 * sin_theta2));

        SVector::<f64, 8>::from_vec(vec![u0, u1, u2, u3, u0_dot, u1_dot, u2_dot, u3_dot])
    }
}

#[cfg(test)]
mod tests {
    use std::{f64::consts::PI, fs::File, io::Write};

    use super::*;

    #[test]
    fn test_kerr() {
        // let incoming: SVector<f64, 8> = SVector::from_vec(vec![
        //     0.,
        //     74.540472912548168,
        //     0.31358562630119541,
        //     2.9957233009687654,
        //     0.050684212954921752,
        //     -0.0029141075809195824,
        //     2.8703288549457436e-05,
        //     0.0021685282189873650,
        // ]);

        let incoming: SVector<f64, 8> = SVector::from_vec(vec![
            0.,
            56.564132427828916,
            0.31358562630119541,
            3.1415926535897931,
            40.720908013320951,
            0.0000000000000000,
            0.0000000000000000,
            1.0000000000000000,
        ]);

        let metric = Metric::new(Tup(0., 0., 0.), 0.999, 2.5);

        let result: SVector<f64, 8> = metric.kerr(incoming);
        let result2: SVector<f64, 8> = metric.derivatives(incoming);
        println!("{:?}", result);
        println!("{:?}", result2);
    }

    #[test]
    fn some_function_test() {
        let camera_position = Tup(0., 0., -10.);
        let num_rays = 40;
        // let ray_range = 10.;
        let h = 0.1;
        let steps = 1000;

        let mut output = File::create("ray_paths.csv").expect("Unable to create file");

        let m = Metric::new(Tup(0., 0., 0.), 0.999, 2.5);

        for i in 0..num_rays {
            let theta = i as f64 * PI / (num_rays as f64 - 1.0);

            // Calculate the direction in the xz-plane
            let x = theta.cos(); // x component (cosine of the angle)
            let z = theta.sin(); // z component (sine of the angle)

            let cam =
                Geodesic::init_cam(camera_position, Tup(0., 0., 1.), Tup(0., 0., 0.), m.clone());

            let geo = Geodesic::init_ray(camera_position, Tup(x, 0., z), &cam);

            let mut path = vec![geo.ray.o];
            let mut current_ray = geo.clone();

            for _ in 0..steps {
                current_ray = geo.m.rk4_kerr(&current_ray, h);
                path.push(current_ray.ray.o);
            }

            // Write path to file
            for point in path {
                writeln!(output, "{},{}", point.0, point.2).expect("Unable to write to file");
            }

            writeln!(output).expect("Unable to write newline");
        }

        println!("Ray paths written to ray_paths.csv");
    }

    #[test]
    fn some_other_function_test() {
        let num_rays = 20;
        let ray_range = 40.;
        let h = 0.1;
        let steps = 9000;

        let mut output = File::create("ray_paths_3d.csv").expect("Unable to create file");

        let m = Metric::new(Tup(0., 0., 0.), -0.999, 2.5);

        for i in 0..num_rays {
            // calculate x from -15 to 15 in equal spaces
            let x = -ray_range + i as f64 * ray_range * 2. / (num_rays as f64 - 1.0);

            let cam = Geodesic::init_cam(
                Tup(x, 0., -40.),
                Tup(0., 0., 1.),
                Tup(0., 0., 0.),
                m.clone(),
            );

            let geo = Geodesic::init_ray(cam.ray.o, Tup(0., 0., 1.), &cam);

            let mut path = vec![geo.ray.o];
            let mut current_ray = geo.clone();

            for _ in 0..steps {
                current_ray = m.rk4_kerr(&current_ray, h);
                path.push(current_ray.ray.o);
            }

            // Write path to file
            for point in path {
                writeln!(output, "{},{},{}", point.0, point.1, point.2)
                    .expect("Unable to write to file");
            }

            writeln!(output).expect("Unable to write newline");
        }

        println!("Ray paths written to ray_paths.csv");
    }
}
