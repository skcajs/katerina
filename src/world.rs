use core::f64;

use super::interval::Metric;
use super::ray::Ray;
use super::sphere::{RflType, Sphere};
use super::tup::Tup;

pub struct World {
    pub spheres: Vec<Sphere>,
}

impl World {
    pub fn new() -> Self {
        World {
            spheres: vec![
                // Scene: radius, position, emission, color, material
                Sphere::new(
                    1e5,
                    Tup(1e5 - 49., -11.2, 81.6 - 25.),
                    Tup::zeros(),
                    Tup(0.75, 0.25, 0.25),
                    RflType::DIFF,
                ), // Left
                Sphere::new(
                    1e5,
                    Tup(-1e5 + 49., -11.2, 81.6 - 25.),
                    Tup::zeros(),
                    Tup(0.25, 0.25, 0.75),
                    RflType::DIFF,
                ), // Right
                Sphere::new(
                    1e5,
                    Tup(0., -11.2, 1e5 - 25.),
                    Tup::zeros(),
                    Tup(0.75, 0.75, 0.75), //Tup(0.25, 0.75, 0.75),
                    RflType::DIFF,
                ), // Back
                Sphere::new(
                    1e5,
                    Tup(0., -11.2, -1e5 + 170. - 25.),
                    Tup::zeros(),
                    Tup::zeros(), //Tup(0.75, 0.75, 0.25),
                    RflType::DIFF,
                ), // Front
                Sphere::new(
                    1e5,
                    Tup(0., 1e5 - 52., 81.6 - 25.),
                    Tup::zeros(),
                    Tup(0.75, 0.75, 0.75),
                    RflType::DIFF,
                ), // Bottom
                Sphere::new(
                    1e5,
                    Tup(0., -1e5 + 29.6, 81.6 - 25.),
                    Tup::zeros(),
                    Tup(0.75, 0.75, 0.75),
                    RflType::DIFF,
                ), // Top
                Sphere::new(
                    16.5,
                    Tup(-23., -35.5, 47.0 - 25.),
                    Tup::zeros(),
                    Tup(1., 1., 1.) * 0.999,
                    RflType::SPEC,
                ), // Mirror
                Sphere::new(
                    16.5,
                    Tup(23., -35.5, 78. - 25.),
                    Tup::zeros(),
                    Tup(1., 1., 1.) * 0.999,
                    RflType::REFR,
                ), // Glass
                Sphere::new(
                    600.,
                    Tup(0., 629.6 - 0.27, 81.6 - 25.),
                    Tup(12., 12., 12.),
                    Tup::zeros(),
                    RflType::DIFF,
                ), // Light
            ],
        }
    }

    pub fn intersect(&self, ray: &Ray, t: &mut f64, id: &mut usize) -> bool {
        *t = f64::INFINITY;
        for i in (0..self.spheres.len()).rev() {
            let d = self.spheres[i].intersect(ray);
            if d != 0.0 && d < *t {
                *t = d;
                *id = i;
            }
        }
        *t < f64::INFINITY
    }

    pub fn trace_geodesic(
        &self,
        ray: &mut Ray,
        t: &mut f64,
        id: &mut usize,
        test_color: &mut bool,
    ) -> bool {
        *t = f64::INFINITY;
        let max_iter = 300.0;
        let mut step_size: f64 = 20.;
        let sigma = 1e-1;

        // Testing
        let metric = Metric::new(5., Tup(0., -11.2, 60.));

        *ray = Ray { o: ray.o, d: ray.d };

        for _ in 0..max_iter as usize {
            let mut hit = false;
            for (i, sphere) in self.spheres.iter().enumerate() {
                let d = sphere.intersect(&ray);
                if d > 0.0 && d < step_size {
                    if d < sigma {
                        *t = d;
                        *id = i;
                        return true;
                    } else {
                        hit = true;
                        step_size *= 0.5;
                        break;
                    }
                }
            }

            if !hit {
                *ray = metric.schwarzschild(ray, step_size);

                // let ray_t = ray.o - metric.s;

                // if metric.rs > 0. && ray_t.1.abs() < 0.1 && ray_t.len() > 10. && ray_t.len() < 30. {
                //     *test_color = true;
                //     return true;
                // }

                // if ray_t.len() <= metric.rs {
                //     return false;
                // }

                step_size = (step_size * 1.05).min(1.0);
            }
        }

        *t = f64::INFINITY;
        false // No intersection found within max_distance
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn new_world() {
        let sphere = Sphere {
            r: 1e5,
            p: Tup(1e5 + 1.0, 40.8, 81.6),
            e: Tup(0., 0., 0.),
            c: Tup(0.75, 0.25, 0.25),
            rfl: RflType::DIFF,
        };
        let world = World::new();
        assert_eq!(world.spheres.len(), 9);
        assert_eq!(world.spheres[0], sphere)
    }

    #[test]
    fn ray_intersects() {
        let world = World::new();
        let ray = Ray {
            o: Tup(0.0, 0.0, -5.0),
            d: Tup(0.0, 0.0, 1.0),
        };
        let mut t = f64::INFINITY;
        let mut id = 0;
        let its = world.intersect(&ray, &mut t, &mut id);
        assert!(its);
    }

    #[test]
    fn ray_does_not_intersect() {
        let world = World::new();
        let ray = Ray {
            o: Tup(0.0, 0.0, -200000.0),
            d: Tup(0.0, 0.0, 0.0),
        };
        let mut t = f64::INFINITY;
        let mut id = 0;
        let its = world.intersect(&ray, &mut t, &mut id);
        assert!(!its);
    }
}
