use super::ray::Ray;
use super::tup::Tup;

#[derive(Debug, PartialEq)]
pub enum RflType {
    DIFF,
    SPEC,
    REFR,
}

#[derive(Debug, PartialEq)]
pub struct Sphere {
    pub r: f64,
    pub p: Tup,
    pub e: Tup,
    pub c: Tup,
    pub rfl: RflType,
}

impl Sphere {
    pub fn new(r: f64, p: Tup, e: Tup, c: Tup, rfl: RflType) -> Self {
        Sphere { r, p, e, c, rfl }
    }

    pub fn intersect(&self, ray: &Ray) -> f64 {
        let eps = 1e-3;
        let op = self.p - ray.o;
        let b = op.dot(ray.d);
        let det = b * b - op.dot(op) + self.r * self.r;
        if det < 0.0 {
            return 0.0;
        }

        let det_sqrt = det.sqrt();
        let mut t = b - det_sqrt;

        if t > eps {
            return t;
        }

        t = b + det_sqrt;

        if t > eps {
            return t;
        }

        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sphere() {
        let r = 1.0;
        let p = Tup(1.0, 2.0, 3.0);
        let e = Tup(1.0, 1.0, 1.0);
        let c = Tup(1.0, 1.0, 1.0);

        let sphere = Sphere {
            r,
            p,
            e,
            c,
            rfl: RflType::DIFF,
        };

        assert_eq!(sphere.r, r);
        assert_eq!(sphere.p, p);
        assert_eq!(sphere.e, e);
        assert_eq!(sphere.c, c);
        assert_eq!(sphere.rfl, RflType::DIFF);
    }

    #[test]
    fn ray_intersects_sphere() {
        let r = 1.0;
        let p = Tup(0.0, 0.0, 0.0);
        let e = Tup(1.0, 1.0, 1.0);
        let c = Tup(1.0, 1.0, 1.0);

        let sphere = Sphere {
            r,
            p,
            e,
            c,
            rfl: RflType::DIFF,
        };

        let ray = Ray {
            o: Tup(0.0, 0.0, -5.0),
            d: Tup(0.0, 0.0, 1.0),
        };

        let xs = sphere.intersect(&ray);
        assert_eq!(xs, 4.0);
    }

    #[test]
    fn ray_misses_sphere() {
        let r = 1.0;
        let p = Tup(0.0, 0.0, 0.0);
        let e = Tup(1.0, 1.0, 1.0);
        let c = Tup(1.0, 1.0, 1.0);

        let sphere = Sphere {
            r,
            p,
            e,
            c,
            rfl: RflType::DIFF,
        };

        let ray = Ray {
            o: Tup(0.0, 2.0, -5.0),
            d: Tup(0.0, 0.0, 1.0),
        };

        let xs = sphere.intersect(&ray);
        assert_eq!(xs, 0.0);
    }
}
