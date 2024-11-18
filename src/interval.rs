use super::ray::Ray;
use super::tup::Tup;

// const G: f64 = 1.;
// const M: f64 = 2.5;
// const RS: f64 = 2. * G * M;
// const RS: f64 = 1.6;

// const S: Tup = Tup(23., -35.5, 78. - 25.);
//122.8

const RS: f64 = 5.;
const S: Tup = Tup(0., -11.5, 60.);

#[allow(dead_code)]
pub fn minkowski(ray: &Ray, h: f64) -> Ray {
    let current_point = ray.o + ray.d * h;

    Ray {
        o: current_point,
        d: (current_point - ray.o).norm(),
    }
}

pub fn schwarzschild(ray: &Ray, h: f64) -> Ray {
    let previous_x = ray.o;
    let previous_p = ray.d;

    let k1x = fx(previous_p, previous_x);
    let k1p = fp(previous_p, previous_x);

    let k2x = fx(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);
    let k2p = fp(previous_p + k1p * 0.5 * h, previous_x + k1x * 0.5 * h);

    let k3x = fx(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);
    let k3p = fp(previous_p + k2p * 0.5 * h, previous_x + k2x * 0.5 * h);

    let k4x = fx(previous_p + k3p * h, previous_x + k3x * h);
    let k4p = fp(previous_p + k3p * h, previous_x + k3x * h);

    let current_point = previous_x + ((k1x + k2x * 2. + k3x * 2. + k4x) * (h / 6.));
    let current_momentum = previous_p + ((k1p + k2p * 2. + k3p * 2. + k4p) * (h / 6.));

    Ray {
        o: current_point,
        d: current_momentum.norm(),
    }
}

pub fn fx(p: Tup, x: Tup) -> Tup {
    let r: f64 = (x - S).len();
    let r_adjusted = r * (1. + (RS / (4. * r))).powi(2);
    let a: f64 = 1. + (RS / (4. * r_adjusted));
    let b: f64 = 1. - (RS / (4. * r_adjusted));
    let fact_x: f64 = (b * b) / a.powi(6);
    p * fact_x
}

pub fn fp(p: Tup, x: Tup) -> Tup {
    let r: f64 = (x - S).len();
    let r_adjusted = r * (1. + (RS / (4. * r))).powi(2);
    let a: f64 = 1. + (RS / (4. * r));
    let b: f64 = 1. - (RS / (4. * r));
    let fact_p1: f64 = (b * b) / a.powi(7);
    let fact_p2: f64 = 1.0 / (b * a);
    x * (-1. / (2. * r_adjusted.powi(3)))
        * (((p.0.powi(2) + p.1.powi(2) + p.2.powi(2)) * fact_p1) + fact_p2)
        * RS
}
