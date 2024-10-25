use super::ray::Ray;
use super::tup::Tup;

const G: f64 = 1.;
const M: f64 = 1.25;
const RS: f64 = 2.*G*M;


pub fn minkowski(ray: &Ray, h: f64) -> Ray {
    let current_point = ray.o + ray.d * h;

    Ray {
        o: current_point,
        d: (current_point - ray.o).norm(),
    }
}

pub fn schwarszchild(ray: &Ray, h: f64) -> Ray {

    let previous_point = ray.o;
    let previous_momentum = ray.d;

    let r: f64 = previous_point.len();

    let a: f64 = 1. + (RS / (4. * r));
    let b: f64 = 1. - (RS / (4. * r));

    let fact_x: f64 = (b * b) / (a * a * a * a * a * a);
    let fact_p1: f64 = (b * b) / (a * a * a * a * a * a * a);
    let fact_p2: f64 = 1.0 / (b * a);

    let k1x = fx(previous_momentum, fact_x);
    let k1p = fp(previous_momentum, previous_point, fact_p1, fact_p2, r);

    let k2x = fx(previous_momentum + k1p * 0.5 * h, fact_x);
    let k2p = fp(
        previous_momentum + k1p * 0.5 * h,
        previous_point + k1x * 0.5 * h,
        fact_p1, fact_p2, r
    );

    let k3x = fx(previous_momentum + k2p * 0.5 * h, fact_x);
    let k3p = fp(
        previous_momentum + k2p * 0.5 * h,
        previous_point + k2x * 0.5 * h,
        fact_p1, fact_p2, r
    );

    let k4x = fx(previous_momentum + k3p * 0.5 * h, fact_x);
    let k4p = fp(previous_momentum + k3p * h, previous_point + k3x * h, fact_p1, fact_p2, r);

    let current_point = previous_point + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
    let current_momentum = previous_momentum + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

    Ray {
        o: current_point,
        d: current_momentum,
    }
}

pub fn fx(p: Tup, fact_x: f64) -> Tup {
    return p * fact_x;
}

pub fn fp(p: Tup, x: Tup, fact_p1: f64, fact_p2: f64, r: f64) -> Tup {
    return x
        * (-1. / (2. * r.powf(3.)))
        * (((p.0.powf(2.) + p.1.powf(2.) + p.2.powf(2.)) * fact_p1) + fact_p2)
        * RS;
}
