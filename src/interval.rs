use num_traits::Pow;

use super::ray::Ray;
use super::tup::Tup;

const G: f64 = 1.;
const M: f64 = 0.;
const R: f64 = 5.;

const RS: f64 = 2.*G*M;
const A: f64 = 1. + (RS / (4. * R));
const B: f64 = 1. - (RS / (4. * R));

const FACT_X: f64 = (B * B) / (A * A * A * A * A * A);
const FACT_P1: f64 = (B * B) / (A * A * A * A * A * A * A);
const FACT_P2: f64 = 1.0 / (B * A);

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

    let k1x = fx(previous_momentum);
    let k1p = fp(previous_momentum, previous_point);

    let k2x = fx(previous_momentum + k1p * 0.5 * h);
    let k2p = fp(
        previous_momentum + k1p * 0.5 * h,
        previous_point + k1x * 0.5 * h,
    );

    let k3x = fx(previous_momentum + k2p * 0.5 * h);
    let k3p = fp(
        previous_momentum + k2p * 0.5 * h,
        previous_point + k2x * 0.5 * h,
    );

    let k4x = fx(previous_momentum + k3p * 0.5 * h);
    let k4p = fp(previous_momentum + k3p * h, previous_point + k3x * h);

    let current_point = previous_point + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
    let current_momentum = previous_momentum + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

    Ray {
        o: current_point,
        d: current_momentum,
    }
}

pub fn fx(p: Tup) -> Tup {
    return p * FACT_X;
}

pub fn fp(p: Tup, x: Tup) -> Tup {
    return x
        * (-1. / (2. * R.pow(3)))
        * (((p.0.pow(2) + p.1.pow(2) + p.2.pow(2)) * FACT_P1) + FACT_P2)
        * RS;
}
