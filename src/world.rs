use core::f64;

use crate::scenes::cornell::cornell_box;

use super::interval::Metric;
use super::ray::Ray;
use super::sphere::Sphere;

pub struct World {
    pub spheres: Vec<Sphere>,
}

impl World {
    pub fn new() -> Self {
        World {
            spheres: cornell_box(),
        }
    }

    #[allow(dead_code)]
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
        accretion_disk: &mut bool,
        metric: &Metric,
    ) -> bool {
        *t = f64::INFINITY;
        let max_iter = 300.0;
        let initial_step_size = 10.;
        let mut step_size: f64 = initial_step_size;
        let sigma = 1e-1;

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
                *ray = metric.rk4(ray, step_size);
                // (*ray, step_size) = metric.rkf45(ray, step_size, 1e-6);

                // println!("{:?}", step_size);

                // let ray_t = ray.o - metric.s;

                // if metric.rs > 0. && ray_t.1.abs() < 0.1 && ray_t.len() > 10. && ray_t.len() < 30. {
                //     *test_color = true;
                //     return true;
                // }

                step_size = (step_size * 1.05).min(initial_step_size * 2.);
            }
        }

        *t = f64::INFINITY;
        false // No intersection found within max_distance
    }
}

#[cfg(test)]
mod tests {

    use crate::{sphere::RflType, tup::Tup};

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
