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
    pub fn rkf45(&self, ray: &Ray, &mut h: &mut f64) -> Ray {
        let previous_x = ray.o - self.s;
        let previous_v = ray.d;

        let mut h_u = h;
        let mut abs_error = 1e10;
        let sigma = 1e-6;

        let incoming: Matrix3<f64> = Matrix3::new(
            0.,
            previous_x.0,
            previous_x.1,
            previous_x.2,
            0.,
            previous_v.0,
            previous_v.1,
            previous_v.2,
            h_u,
        );

        let mut y: Matrix3<f64> = Matrix3::zeros();

        let max_iterations = 2;
        let mut iterations = 0;
        while abs_error > sigma && iterations < max_iterations {
            iterations += 1;
            let k1: Matrix3<f64> = self.derivatives(incoming);
            let k2: Matrix3<f64> = self.derivatives(incoming + k1 * (h_u / 4.));
            let k3: Matrix3<f64> =
                self.derivatives(incoming + k1 * h_u * (3. / 32.) + k2 * h_u * (9. / 32.));
            let k4: Matrix3<f64> = self.derivatives(
                incoming + k1 * h_u * (1932. / 2197.) - k2 * h_u * (7200. / 2197.)
                    + k3 * h_u * (7296. / 2197.),
            );
            let k5: Matrix3<f64> = self.derivatives(
                incoming + k1 * h_u * (439. / 216.) - k2 * h_u * 8. + k3 * h_u * (3680. / 513.)
                    - k4 * h_u * (845. / 4104.),
            );
            let k6: Matrix3<f64> = self.derivatives(
                incoming - k1 * h_u * (8. / 27.) + k2 * h_u * 2. - k3 * h_u * (3544. / 2565.)
                    + k4 * h_u * (1859. / 4104.)
                    - k5 * h_u * (11. / 40.),
            );

            y = incoming
                + (k1 * h_u * (16. / 135.)
                    + k3 * h_u * (6656. / 12825.)
                    + k4 * h_u * (28561. / 56430.)
                    - k5 * h_u * (9. / 50.)
                    + k6 * h_u * (2. / 55.));

            let error: Matrix3<f64> =
                k1 * h_u * (25. / 216.) + k3 * h_u * (1408. / 2565.) + k4 * h_u * (2197. / 4104.)
                    - k5 * h_u * (1. / 5.);

            abs_error = f64::sqrt(
                error.m12.powi(2)
                    + error.m13.powi(2)
                    + error.m21.powi(2)
                    + error.m23.powi(2)
                    + error.m31.powi(2)
                    + error.m32.powi(2),
            );

            let new_h = h_u * (sigma / abs_error).powf(0.2);

            if abs_error > sigma {
                h_u = new_h;
            }
        }

        let current_p = Tup(y.m12, y.m13, y.m21);
        let current_v = Tup(y.m23, y.m31, y.m32).norm();

        Ray {
            o: current_p + self.s,
            d: current_v,
        }
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

    fn derivatives(&self, y: Matrix3<f64>) -> Matrix3<f64> {
        let x = Tup(y.m12, y.m13, y.m21);
        let r: f64 = x.len();
        let p: Tup = Tup(y.m23, y.m31, y.m32);
        if r < self.rs / 4. {
            return Matrix3::new(0., y.m12, y.m13, y.m21, 0., y.m23, y.m31, y.m32, 0.);
        }
        let r_adjusted = r * (1. + (self.rs / (4. * r))).powi(2);
        let a: f64 = 1. + (self.rs / (4. * r));
        let b: f64 = 1. - (self.rs / (4. * r));
        let fact_x: f64 = (b * b) / a.powi(6);
        let fact_p1: f64 = (b * b) / a.powi(7);
        let fact_p2: f64 = 1.0 / (b * a);
        let p_new = p * fact_x;
        let x_new = x
            * (-1. / (2. * r_adjusted.powi(3)))
            * (((p.0.powi(2) + p.1.powi(2) + p.2.powi(2)) * fact_p1) + fact_p2)
            * self.rs;

        Matrix3::new(
            0., x_new.0, x_new.1, x_new.2, 0., p_new.0, p_new.1, p_new.2, 0.,
        )
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
