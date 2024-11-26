use recursive::recursive;
use std::f64::consts::PI;

use crate::{
    geodesic::Geodesic,
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
    ray: Geodesic,
    depth: i32,
    sampler: &mut Sampler,
    int_type: IntegrationType,
) -> Tup {
    match int_type {
        IntegrationType::Iterative => radiance_iter(world, ray, depth, sampler),
        IntegrationType::Recursive => radiance(world, ray, depth, sampler, 1.),
    }
}

pub fn radiance_iter(
    world: &World,
    initial_ray: Geodesic,
    depth: i32,
    sampler: &mut Sampler,
) -> Tup {
    // let metric = Metric::new(1.6, Tup(23., -35.5, 78. - 25.));
    // let metric = Metric::new(5., Tup(-1., -13.2, 60.));
    let mut result = Tup::zeros();
    let throughput = Tup::ones();
    let mut stack = Vec::new(); // Stack to handle rays iteratively
    stack.push((initial_ray, throughput, depth));

    while let Some((mut geo, mut throughput, mut depth)) = stack.pop() {
        loop {
            let mut t = f64::INFINITY;
            let mut id: usize = 0;
            let mut accretion_disk = false;

            if !world.trace_geodesic(&mut geo, &mut t, &mut id, &mut accretion_disk) {
                break;
            }

            // if accretion_disk {
            //     result += Tup(0.8, 0., 0.8) * throughput;
            //     break;
            // }

            let obj: &Sphere = &world.spheres[id];
            let x = geo.ray.o + (geo.ray.d * t);
            let n = (x - obj.p).norm();
            let n1 = if n.dot(geo.ray.d) < 0.0 { n } else { n * -1.0 };

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
                    geo = Geodesic::ray(x, d, geo.m);
                }
                RflType::SPEC => {
                    geo = Geodesic::ray(x, geo.ray.d - n * 2. * n.dot(geo.ray.d), geo.m);
                }
                RflType::REFR => {
                    let into = n.dot(n1) > 0.;
                    let nc: f64 = 1.;
                    let nt: f64 = 1.5;
                    let nnt = if into { nc / nt } else { nt / nc };
                    let ddn = geo.ray.d.dot(n1);
                    let cos2t = 1. - nnt * nnt * (1. - ddn * ddn);

                    if cos2t < 0. {
                        geo = Geodesic::ray(x, geo.ray.d - n * 2. * n.dot(geo.ray.d), geo.m);
                        continue;
                    }

                    let tdir = (geo.ray.d * nnt
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
                            geo = Geodesic::ray(x, geo.ray.d - n * 2. * n.dot(geo.ray.d), geo.m);
                            throughput = throughput * rp;
                        } else {
                            geo.ray = Ray { o: x, d: tdir };
                            throughput = throughput * tp;
                        }
                    } else {
                        let reflected_dir = geo.ray.d - n * 2. * n.dot(geo.ray.d);
                        let reflected_throughput = throughput * re;
                        let refracted_throughput = throughput * tr;

                        stack.push((
                            Geodesic::ray(x, reflected_dir, geo.m.clone()),
                            reflected_throughput,
                            depth,
                        ));
                        stack.push((
                            Geodesic::ray(x, tdir, geo.m.clone()),
                            refracted_throughput,
                            depth,
                        ));
                        break;
                    }
                }
            }
        }
    }
    result
}

#[recursive]
pub fn radiance(
    world: &World,
    mut geo: Geodesic,
    mut depth: i32,
    mut sampler: &mut Sampler,
    em: f64,
) -> Tup {
    // let metric = Metric::new(1.6, Tup(23., -35.5, 78. - 25.));
    // let metric = Metric::new(5., Tup(-1., -13.2, 60.));
    let mut t = f64::INFINITY;
    let mut id: usize = 0;
    if !world.trace_geodesic(&mut geo, &mut t, &mut id, &mut false) {
        return Tup::zeros();
    }
    let obj: &Sphere = &world.spheres[id];
    let x = geo.ray.o + (geo.ray.d * t); // I think this should be the current point from the trace_geodesic.
    let n = (x - obj.p).norm();
    let n1 = if n.dot(geo.ray.d) < 0.0 { n } else { n * -1.0 };

    let mut f = obj.c;
    let p = f.0.max(f.1.max(f.2));
    depth += 1;
    if depth > 5 || !p.is_finite() {
        if sampler.next() < p {
            f = f * (1.0 / p);
        } else {
            return obj.e * em;
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

            // Loop over any lights
            let mut e = Tup::zeros();
            for (i, sphere) in world.spheres.iter().enumerate() {
                if !sphere.is_emitter() {
                    continue; // Skip non-emissive objects
                }
                // let adjusted_sphere = metric.transform_point(sphere.p, x);
                // Calculate direction to light
                let sw = (sphere.p - x).norm();
                let su = if sw.0.abs() > 0.1 {
                    Tup(0., 1., 0.).cross(sw).norm()
                } else {
                    Tup(1., 0., 0.).cross(sw).norm()
                };
                let sv = sw.cross(su);

                let cos_a_max =
                    (1. - sphere.r * sphere.r / (x - sphere.p).dot(x - sphere.p)).sqrt();

                let eps1 = sampler.next();
                let eps2 = sampler.next();
                let cos_a = 1. - eps1 + eps1 * cos_a_max;
                let sin_a = (1. - cos_a * cos_a).sqrt();
                let phi = 2. * PI * eps2;
                let l =
                    (su * f64::cos(phi) * sin_a + sv * f64::sin(phi) * sin_a + sw * cos_a).norm();
                if world.trace_geodesic(
                    &mut Geodesic::ray(x, l, geo.m.clone()),
                    &mut t,
                    &mut id,
                    &mut false,
                ) && id == i
                {
                    let omega = 2. * PI * (1. - cos_a_max);
                    e += f * sphere.e * f64::max(0., l.dot(n1)) * omega * (1. / PI);
                }
            }

            return obj.e * em
                + e
                + (f * radiance(
                    world,
                    Geodesic::ray(x, d, geo.m.clone()),
                    depth,
                    &mut sampler,
                    0.,
                ));
        }
        RflType::SPEC => {
            return obj.e
                + f * radiance(
                    world,
                    Geodesic::ray(x, geo.ray.d - n * 2. * n.dot(geo.ray.d), geo.m),
                    depth,
                    &mut sampler,
                    1.,
                );
        }
        RflType::REFR => {
            let rfl_ray = Geodesic::ray(x, geo.ray.d - n * 2. * n.dot(geo.ray.d), geo.m.clone());
            let into = n.dot(n1) > 0.;
            let nc: f64 = 1.;
            let nt: f64 = 1.5;
            let nnt = if into { nc / nt } else { nt / nc };
            let ddn = geo.ray.d.dot(n1);
            let cos2t = 1. - nnt * nnt * (1. - ddn * ddn);
            if cos2t < 0. {
                return obj.e + f * radiance(world, rfl_ray, depth, &mut sampler, 1.);
            }
            let tdir = (geo.ray.d * nnt
                - n * if into { 1. } else { -1. } * (ddn * nnt + cos2t.sqrt()))
            .norm();
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
                        radiance(world, rfl_ray, depth, &mut sampler, 1.) * rp
                    } else {
                        radiance(
                            world,
                            Geodesic::ray(x, tdir, geo.m.clone()),
                            depth,
                            &mut sampler,
                            1.,
                        ) * tp
                    }
                } else {
                    radiance(world, rfl_ray, depth, &mut sampler, 1.) * re
                        + radiance(
                            world,
                            Geodesic::ray(x, tdir, geo.m.clone()),
                            depth,
                            &mut sampler,
                            1.,
                        ) * tr
                })
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::metric::Metric;

    use super::*;
    #[test]
    fn ray_intesects_empty_world() {
        let geo = Geodesic::ray(
            Tup(0., 0., 0.),
            Tup(0., 0., -1.),
            Metric {
                s: Tup(0., 0., 0.),
                a: -0.999,
                rs: 0.0,
            },
        );
        let world = World { spheres: vec![] };
        let mut sampler = Sampler::new();

        let result = radiance(&world, geo, 0, &mut sampler, 1.);
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
        let geo = Geodesic::ray(
            Tup(0., 0., 0.),
            Tup(0., 0., -1.),
            Metric {
                s: Tup(0., 0., 0.),
                a: -0.999,
                rs: 0.0,
            },
        );
        let mut sampler = Sampler::new();

        let result = radiance(&world, geo, 0, &mut sampler, 1.);
        assert_eq!(result, Tup(1., 0., 0.));
    }
}
