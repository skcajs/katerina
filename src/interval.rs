extern crate nalgebra as na;

use na::Matrix3;

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
    pub fn rkf45(&self, ray: &Ray, h: f64, tol: f64) -> (Ray, f64) {
        let mut h = h;
        let previous_x = ray.o - self.s;
        let previous_p = ray.d;

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
                return (
                    Ray {
                        o: x5 + self.s,
                        d: p5.norm(),
                    },
                    h,
                );
            }

            // Reduce step size for better accuracy
            h *= 0.9 * (tol / error).powf(0.2);
        }
    }

    #[allow(dead_code)]
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
