use super::ray::Ray;
use super::tup::Tup;

#[derive(Debug, Clone, Copy)]
pub struct Metric {
    pub rs: f64,
    pub s: Tup,
}

impl Metric {
    pub fn new(rs: f64, s: Tup) -> Self {
        Metric { rs, s }
    }

    #[allow(dead_code)]
    pub fn minkowski(&self, ray: &Ray, h: f64) -> Ray {
        let current_point = ray.o + ray.d * h;

        Ray {
            o: current_point,
            d: (current_point - ray.o).norm(),
        }
    }

    pub fn schwarzschild(&self, ray: &Ray, h: f64) -> Ray {
        let previous_x = ray.o - self.s;
        let previous_p = ray.d;

        let k1x = self.fx(previous_p, previous_x);
        let k1p = self.fp(previous_p, previous_x);

        let k2x = self.fx(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);
        let k2p = self.fp(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);

        let k3x = self.fx(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);
        let k3p = self.fp(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);

        let k4x = self.fx(previous_p + k3p * h, previous_x + k3x * h);
        let k4p = self.fp(previous_p + k3p * h, previous_x + k3x * h);

        let current_point = previous_x + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
        let current_momentum = previous_p + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

        Ray {
            o: current_point + self.s,
            d: current_momentum.norm(),
        }
    }

    fn fx(&self, p: Tup, x: Tup) -> Tup {
        let r: f64 = x.len();
        // let r_adjusted = r * (1. + (RS / (4. * r))).powi(2);
        let a: f64 = 1. + (self.rs / (4. * r));
        let b: f64 = 1. - (self.rs / (4. * r));
        let fact_x: f64 = (b * b) / a.powi(6);
        p * fact_x
    }

    fn fp(&self, p: Tup, x: Tup) -> Tup {
        let r: f64 = x.len();
        // let r_adjusted = r * (1. + (RS / (4. * r))).powi(2);
        let a: f64 = 1. + (self.rs / (4. * r));
        let b: f64 = 1. - (self.rs / (4. * r));
        let fact_p1: f64 = (b * b) / a.powi(7);
        let fact_p2: f64 = 1.0 / (b * a);
        x * (-1. / (2. * r.powi(3)))
            * (((p.0.powi(2) + p.1.powi(2) + p.2.powi(2)) * fact_p1) + fact_p2)
            * self.rs
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

        let metric = Metric::new(2.5, Tup(0., -11.2, 65.));

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
                current_ray = metric.schwarzschild(&current_ray, h);
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

        let metric = Metric::new(2.5, Tup(0., 0., 0.));

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
