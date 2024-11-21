use recursive::recursive;
use std::f64::consts::PI;

use crate::{
    ray::Ray,
    sampler::Sampler,
    sphere::{RflType, Sphere},
    tup::Tup,
    world::World,
};

#[derive(Default)]
pub enum IntegrationType {
    #[default]
    Iterative,
    #[allow(dead_code)]
    Recursive,
}

pub fn integrate(
    world: &World,
    ray: Ray,
    depth: i32,
    sampler: &mut Sampler,
    int_type: IntegrationType,
) -> Tup {
    match int_type {
        IntegrationType::Iterative => radiance_iter(world, ray, depth, sampler),
        IntegrationType::Recursive => radiance(world, ray, depth, sampler),
    }
}

pub fn radiance_iter(world: &World, initial_ray: Ray, depth: i32, sampler: &mut Sampler) -> Tup {
    let mut result = Tup::zeros();
    let throughput = Tup::ones();
    let mut stack = Vec::new(); // Stack to handle rays iteratively
    stack.push((initial_ray, throughput, depth));

    while let Some((mut ray, mut throughput, mut depth)) = stack.pop() {
        loop {
            let mut t = f64::INFINITY;
            let mut id: usize = 0;
            let mut test_color = false;

            if !world.trace_geodesic(&mut ray, &mut t, &mut id, &mut test_color) {
                break;
            }

            // if test_color {
            //     result += Tup(0.8, 0., 0.8) * throughput;
            //     break;
            // }

            let obj: &Sphere = &world.spheres[id];
            let x = ray.o + (ray.d * t);
            let n = (x - obj.p).norm();
            let n1 = if n.dot(ray.d) < 0.0 { n } else { n * -1.0 };

            let mut f = obj.c;
            let p = f.0.max(f.1.max(f.2));
            depth += 1;

            if depth > 5 {
                if sampler.next() < p {
                    f = f * (1.0 / p);
                } else {
                    result += throughput * obj.e;
                    break;
                }
            }

            result += throughput * obj.e;
            throughput = throughput * f;

            match obj.rfl {
                RflType::DIFF => {
                    let r1 = 2. * PI * sampler.next();
                    let r2: f64 = sampler.next();
                    let r2s = r2.sqrt();
                    let w: Tup = n1;
                    let u: Tup = if w.0.abs() > 0.1 {
                        Tup(0., 1., 0.).cross(w).norm()
                    } else {
                        Tup(1., 0., 0.).cross(w).norm()
                    };
                    let v = w.cross(u);
                    let d: Tup =
                        (u * f64::cos(r1) * r2s + v * f64::sin(r1) * r2s + w * ((1. - r2).sqrt()))
                            .norm();

                    ray = Ray { o: x, d };
                }
                RflType::SPEC => {
                    ray = Ray {
                        o: x,
                        d: ray.d - n * 2. * n.dot(ray.d),
                    };
                }
                RflType::REFR => {
                    let into = n.dot(n1) > 0.;
                    let nc: f64 = 1.;
                    let nt: f64 = 1.5;
                    let nnt = if into { nc / nt } else { nt / nc };
                    let ddn = ray.d.dot(n1);
                    let cos2t = 1. - nnt * nnt * (1. - ddn * ddn);

                    if cos2t < 0. {
                        ray = Ray {
                            o: x,
                            d: ray.d - n * 2. * n.dot(ray.d),
                        };
                        continue;
                    }

                    let tdir = (ray.d * nnt
                        - n * if into { 1. } else { -1. } * (ddn * nnt + cos2t.sqrt()))
                    .norm();

                    let a = nt - nc;
                    let b = nt + nc;
                    let r0 = (a * a) / (b * b);
                    let c = 1. - if into { -ddn } else { tdir.dot(n) };
                    let re = r0 + (1. - r0) * c * c * c * c * c;
                    let tr = 1. - re;
                    let p = 0.25 + 0.5 * re;
                    let rp = re / p;
                    let tp = tr / (1. - p);

                    if depth > 2 {
                        if sampler.next() < p {
                            ray = Ray {
                                o: x,
                                d: ray.d - n * 2. * n.dot(ray.d),
                            };
                            throughput = throughput * rp;
                        } else {
                            ray = Ray { o: x, d: tdir };
                            throughput = throughput * tp;
                        }
                    } else {
                        let reflected_dir = ray.d - n * 2. * n.dot(ray.d);
                        let reflected_throughput = throughput * re;
                        let refracted_throughput = throughput * tr;

                        stack.push((
                            Ray {
                                o: x,
                                d: reflected_dir,
                            },
                            reflected_throughput,
                            depth,
                        ));
                        stack.push((Ray { o: x, d: tdir }, refracted_throughput, depth));
                        break;
                    }
                }
            }
        }
    }
    result
}

#[recursive]
pub fn radiance(world: &World, mut ray: Ray, mut depth: i32, mut sampler: &mut Sampler) -> Tup {
    let mut t = f64::INFINITY;
    let mut id: usize = 0;
    let mut test_color = false;
    if !world.trace_geodesic(&mut ray, &mut t, &mut id, &mut test_color) {
        return Tup(0., 0., 0.);
    }
    if test_color {
        return Tup(0.8, 0., 0.8);
    }
    let obj: &Sphere = &world.spheres[id];
    let x = ray.o + (ray.d * t); // I think this should be the current point from the trace_geodesic.
    let n = (x - obj.p).norm();
    let n1 = if n.dot(ray.d) < 0.0 { n } else { n * -1.0 };

    let mut f = obj.c;
    let p = f.0.max(f.1.max(f.2));
    depth += 1;
    if depth > 5 {
        if sampler.next() < p {
            f = f * (1.0 / p);
        } else {
            return obj.e;
        }
    }

    match obj.rfl {
        RflType::DIFF => {
            let r1 = 2. * PI * sampler.next();
            let r2: f64 = sampler.next();
            let r2s = r2.sqrt();
            let w: Tup = n1;
            let u: Tup = if w.0.abs() > 0.1 {
                Tup(0., 1., 0.).cross(w).norm()
            } else {
                Tup(1., 0., 0.).cross(w).norm()
            };

            let v = w.cross(u);
            let d: Tup =
                (u * f64::cos(r1) * r2s + v * f64::sin(r1) * r2s + w * ((1. - r2).sqrt())).norm();
            return obj.e + f * radiance(world, Ray { o: x, d }, depth, &mut sampler);
        }
        RflType::SPEC => {
            return obj.e
                + f * radiance(
                    world,
                    Ray {
                        o: x,
                        d: ray.d - n * 2. * n.dot(ray.d),
                    },
                    depth,
                    &mut sampler,
                );
        }
        RflType::REFR => {
            let rfl_ray = Ray {
                o: x,
                d: ray.d - n * 2. * n.dot(ray.d),
            };
            let into = n.dot(n1) > 0.;
            let nc: f64 = 1.;
            let nt: f64 = 1.5;
            let nnt = if into { nc / nt } else { nt / nc };
            let ddn = ray.d.dot(n1);
            let cos2t = 1. - nnt * nnt * (1. - ddn * ddn);
            if cos2t < 0. {
                return obj.e + f * radiance(world, rfl_ray, depth, &mut sampler);
            }
            let tdir =
                (ray.d * nnt - n * if into { 1. } else { -1. } * (ddn * nnt + cos2t.sqrt())).norm();
            let a = nt - nc;
            let b = nt + nc;
            let r0 = a * a / (b * b);
            let c = 1. - if into { -ddn } else { tdir.dot(n) };
            let re = r0 + (1. - r0) * c * c * c * c * c;
            let tr = 1. - re;
            let p = 0.25 + 0.5 * re;
            let rp = re / p;
            let tp = tr / (1. - p);

            obj.e
                + f * (if depth > 2 {
                    if sampler.next() < p {
                        radiance(world, rfl_ray, depth, &mut sampler) * rp
                    } else {
                        radiance(world, Ray { o: x, d: tdir }, depth, &mut sampler) * tp
                    }
                } else {
                    radiance(world, rfl_ray, depth, &mut sampler) * re
                        + radiance(world, Ray { o: x, d: tdir }, depth, &mut sampler) * tr
                })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ray_intesects_empty_world() {
        let ray = Ray {
            o: Tup(0., 0., 0.),  // Origin
            d: Tup(0., 0., -1.), // Direction pointing away from any spheres
        };
        let world = World { spheres: vec![] };
        let mut sampler = Sampler::new();

        let result = radiance(&world, ray, 0, &mut sampler);
        assert_eq!(result, Tup(0., 0., 0.));
    }

    #[test]
    fn ray_intesects_single_sphere_world() {
        let sphere = Sphere::new(
            1.0,
            Tup(0., 0., -5.),
            Tup(1., 0., 0.),
            Tup(0., 0., 0.),
            RflType::DIFF,
        );
        let world = World {
            spheres: vec![sphere],
        };
        let ray = Ray {
            o: Tup(0., 0., 0.),  // Origin
            d: Tup(0., 0., -1.), // Direction pointing away from any spheres
        };
        let mut sampler = Sampler::new();

        let result = radiance(&world, ray, 0, &mut sampler);
        assert_eq!(result, Tup(1., 0., 0.));
    }
}
