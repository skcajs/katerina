extern crate nalgebra as na;

use na::Matrix3;

use crate::geodesic::Geodesic;

use super::ray::Ray;
use super::tup::Tup;

#[derive(Debug, Clone, Copy)]
pub struct Metric {
    pub rs: f64,
    pub s: Tup,
}

impl Metric {
    pub fn new(rs: f64, s: Tup) -> Self {
        Metric { rs: rs * 0.6, s } // applying a scale factor to the Schwarzschild radius
    }

    #[allow(dead_code)]
    pub fn minkowski(&self, ray: &Ray, h: f64) -> Ray {
        let current_point = ray.o + ray.d * h;

        Ray {
            o: current_point,
            d: (current_point - ray.o).norm(),
        }
    }

    #[allow(dead_code)]
    pub fn rkf45(&self, incoming: &Geodesic, h: &mut f64) -> Geodesic {
        let mut abs_err = 1e10;
        let tol = 1e-5;
        let mut new_h = 0.0;

        let incoming_vals: Matrix3<f64> = Matrix3::new(
            incoming.t,
            incoming.r,
            incoming.theta,
            incoming.phi,
            incoming.u0,
            incoming.u1,
            incoming.u2,
            incoming.u3,
            incoming.h_step,
        );

        let a = incoming.a;

        let mut h_u = incoming_vals[8];

        let mut new_y: Matrix3<f64> = Matrix3::zeros();

        while abs_err > tol {
            let k1: Matrix3<f64> = self.kerr(incoming_vals, a);
            let k2: Matrix3<f64> = self.kerr(incoming_vals + (1.0 / 4.0) * k1 * h_u, a);
            let k3: Matrix3<f64> = self.kerr(
                incoming_vals + (3.0 / 32.0) * k1 * h_u + (9.0 / 32.0) * k2 * h_u,
                a,
            );
            let k4: Matrix3<f64> = self.kerr(
                incoming_vals + (1932.0 / 2197.0) * k1 * h_u - (7200.0 / 2197.0) * k2 * h_u
                    + (7296.0 / 2197.0) * k3 * h_u,
                a,
            );
            let k5: Matrix3<f64> = self.kerr(
                incoming_vals + (439.0 / 216.0) * k1 * h_u - 8.0 * k2 * h_u
                    + (3680.0 / 513.0) * k3 * h_u
                    - (845.0 / 4104.0) * k4 * h_u,
                a,
            );
            let k6: Matrix3<f64> = self.kerr(
                incoming_vals - (8.0 / 27.0) * k1 * h_u + 2.0 * k2 * h_u
                    - (3544.0 / 2565.0) * k3 * h_u
                    + (1859.0 / 4104.0) * k4 * h_u
                    - (11.0 / 40.0) * k5 * h_u,
                a,
            );

            new_y = incoming_vals
                + (16.0 / 135.0) * k1 * h_u
                + (6656.0 / 12825.0) * k3 * h_u
                + (28561.0 / 56430.0) * k4 * h_u
                - (9.0 / 50.0) * k5 * h_u
                + (2.0 / 55.0) * k6 * h_u;
            let error: Matrix3<f64> = (-1.0 / 360.0) * k1 * h_u
                + (128.0 / 4275.0) * k3 * h_u
                + (2197.0 / 75240.0) * k4 * h_u
                - (1.0 / 50.0) * k5 * h_u
                - (2.0 / 55.0) * k6 * h_u;
            abs_err = f64::sqrt(
                (error[0] * error[0])
                    + (error[1] * error[1])
                    + (error[2] * error[2])
                    + (error[3] * error[3])
                    + (error[4] * error[4])
                    + (error[5] * error[5])
                    + (error[6] * error[6])
                    + (error[7] * error[7])
                    + (error[8] * error[8]),
            );

            new_h = 0.9 * h_u * (tol / abs_err).powf(1.0 / 5.0);

            if abs_err > tol {
                h_u = new_h;
            }
        }

        let geodesic = Geodesic {
            t: new_y[0],
            r: new_y[1],
            theta: new_y[2],
            phi: new_y[3],
            u0: new_y[4],
            u1: new_y[5],
            u2: new_y[6],
            u3: new_y[7],
            h_step: new_h,
            dx: incoming.dx,
            a,
        };
        return geodesic;
    }

    pub fn rk4(&self, ray: &Ray, h: f64) -> Ray {
        // println!("h: {}", h);
        let previous_x = ray.o - self.s;
        let previous_p = ray.d;

        let (k1x, k1p) = self.schwarzschild(previous_p, previous_x);
        let (k2x, k2p) = self.schwarzschild(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);
        let (k3x, k3p) = self.schwarzschild(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);
        let (k4x, k4p) = self.schwarzschild(previous_p + k3p * h, previous_x + k3x * h);

        let current_point = previous_x + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
        let current_momentum = previous_p + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

        Ray {
            o: current_point + self.s,
            d: current_momentum.norm(),
        }
    }

    #[allow(dead_code)]
    pub fn euler(&self, ray: &Ray, h: f64) -> Ray {
        let pos = ray.o;
        let dir = ray.d.norm();

        let r: f64 = (pos - self.s).len();

        if r < 1. {
            return Ray { o: pos, d: dir };
        }

        let h2 = pos.cross(dir).len().powi(2);

        let new_pos = pos + dir * h;
        let new_dir = dir + (pos * -1.5 * h2 * (1. / r.powi(2).powf(2.5) * h));

        Ray {
            o: new_pos,
            d: new_dir,
        }
    }

    fn kerr(&self, incoming_vals: Matrix3<f64>, a: f64) -> Matrix3<f64> {
        let r = incoming_vals[1];
        let theta = incoming_vals[2];
        let u0 = incoming_vals[4];
        let u1 = incoming_vals[5];
        let u2 = incoming_vals[6];
        let u3 = incoming_vals[7];
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

        Matrix3::new(u0, u3, u2_dot, u1, u0_dot, u3_dot, u2, u1_dot, 0.0)
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::{f64::consts::PI, fs::File, io::Write};

    #[test]
    fn some_function_test() {
        let camera_position = Tup(0., 0., -10.);
        let num_rays = 100;
        // let ray_range = 10.;
        let h = 5.;
        let steps = 200;

        let mut output = File::create("ray_paths.csv").expect("Unable to create file");

        let metric = Metric::new(5., Tup(0., 0., 0.));

        for i in 0..num_rays {
            let theta = i as f64 * PI / (num_rays as f64 - 1.0);

            // Calculate the direction in the xz-plane
            let x = theta.cos(); // x component (cosine of the angle)
            let z = theta.sin(); // z component (sine of the angle)

            let ray = Ray {
                o: camera_position,
                d: Tup(x, 0., z),
            };

            let mut path = vec![ray.o];
            let mut current_ray = ray;

            for _ in 0..steps {
                current_ray = metric.rk4(&current_ray, h);
                path.push(current_ray.o);
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
        let num_rays = 40;
        let ray_range = 40.;
        let h = 0.5;
        let steps = 1000;

        let mut output = File::create("ray_paths.csv").expect("Unable to create file");

        let metric = Metric::new(5., Tup(0., 0., 0.));

        for i in 0..num_rays {
            // calculate x from -15 to 15 in equal spaces
            let x = -ray_range + i as f64 * ray_range * 2. / (num_rays as f64 - 1.0);

            let ray = Ray {
                o: Tup(x, 0., -40.),
                d: Tup(0., 0., 1.),
            };

            let mut path = vec![ray.o];
            let mut current_ray = ray;

            for _ in 0..steps {
                current_ray = metric.euler(&current_ray, h);
                path.push(current_ray.o);
            }

            // Write path to file
            for point in path {
                writeln!(output, "{},{}", point.0, point.2).expect("Unable to write to file");
            }

            writeln!(output).expect("Unable to write newline");
        }

        println!("Ray paths written to ray_paths.csv");
    }
}
