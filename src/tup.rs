use std::ops;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Tup(pub f64, pub f64, pub f64);

impl Tup {
    pub fn zeros() -> Self {
        Tup(0., 0., 0.)
    }

    pub fn ones() -> Self {
        Tup(1., 1., 1.)
    }

    pub fn len(&self) -> f64 {
        f64::sqrt(self.0.powf(2.) + (self.1).powf(2.) + self.2.powf(2.))
    }

    pub fn norm(self) -> Self {
        self * (1.0 / (self.0 * self.0 + self.1 * self.1 + self.2 * self.2).sqrt())
    }

    pub fn dot(self, rhs: Tup) -> f64 {
        self.0 * rhs.0 + self.1 * rhs.1 + self.2 * rhs.2
    }

    pub fn cross(self, rhs: Tup) -> Tup {
        Tup(
            self.1 * rhs.2 - self.2 * rhs.1,
            self.2 * rhs.0 - self.0 * rhs.2,
            self.0 * rhs.1 - self.1 * rhs.0,
        )
    }

    #[allow(dead_code)]
    pub fn spherical_to_cartesian(self) -> Self {
        // r = self.0, theta = self.1, phi = self.2
        Tup(
            self.0 * f64::sin(self.1) * f64::cos(self.2),
            self.0 * f64::sin(self.1),
            self.0 * f64::cos(self.1) * f64::sin(self.2),
        )
    }

    #[allow(dead_code)]
    pub fn cartesian_to_spherical(self) -> Self {
        // x = self.0, y = self.1, z = self.2
        let r = self.len();
        let theta = f64::acos(self.1 / r);
        let phi = f64::atan2(self.2, self.0);
        Tup(r, theta, phi)
    }

    #[allow(dead_code)]
    pub fn cartesian_to_boyer_lindquist(self, a: f64) -> Self {
        let Tup(x, y, z) = self;
        let w = (x * x + y * y + z * z) - (a * a);
        let r = f64::sqrt(0.5 * (w + f64::sqrt((w * w) + (4. * (a * a) * (y * y)))));
        let theta = f64::acos(y / r);
        let phi = f64::atan2(z, x);
        Tup(r, theta, phi)
    }

    #[allow(dead_code)]
    pub fn boyer_lindquist_to_cartesian(self, a: f64) -> Self {
        let Tup(r, theta, phi) = self;
        let sqrt_term = (r * r + a * a).sqrt();
        let x = sqrt_term * theta.sin() * phi.cos();
        let y = r * theta.cos();
        let z = sqrt_term * theta.sin() * phi.sin();
        Tup(x, y, z)
    }
}

impl ops::Add<Tup> for Tup {
    type Output = Tup;

    fn add(self, rhs: Tup) -> Self::Output {
        Tup(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl ops::Sub<Tup> for Tup {
    type Output = Tup;

    fn sub(self, rhs: Tup) -> Self::Output {
        Tup(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2)
    }
}

impl ops::Mul<Tup> for Tup {
    type Output = Tup;

    fn mul(self, rhs: Tup) -> Self::Output {
        Tup(self.0 * rhs.0, self.1 * rhs.1, self.2 * rhs.2)
    }
}

impl ops::Mul<f64> for Tup {
    type Output = Tup;

    fn mul(self, rhs: f64) -> Self::Output {
        Tup(self.0 * rhs, self.1 * rhs, self.2 * rhs)
    }
}

impl ops::AddAssign<Tup> for Tup {
    fn add_assign(&mut self, rhs: Tup) {
        *self = Self(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_vec() {
        let v = Tup(1.0, 2.0, 3.0);

        assert_eq!(v.0, 1.0);
        assert_eq!(v.1, 2.0);
        assert_eq!(v.2, 3.0);
    }

    #[test]
    fn vec_add() {
        let v1 = Tup(1.0, 2.0, 3.0);
        let v2 = Tup(2.0, 3.0, 4.0);
        let v3 = v1 + v2;
        assert_eq!(v3.0, 3.0);
        assert_eq!(v3.1, 5.0);
        assert_eq!(v3.2, 7.0);
    }

    #[test]
    fn vec_sub() {
        let v1 = Tup(3.0, 2.0, 1.0);
        let v2 = Tup(1.0, 2.0, 3.0);
        let v3 = v1 - v2;
        assert_eq!(v3.0, 2.0);
        assert_eq!(v3.1, 0.0);
        assert_eq!(v3.2, -2.0);
    }

    #[test]
    fn vec_mul_vec() {
        let v1 = Tup(3.0, 2.0, 1.0);
        let v2 = Tup(1.0, 2.0, 3.0);
        let v3 = v1 * v2;
        assert_eq!(v3.0, 3.0);
        assert_eq!(v3.1, 4.0);
        assert_eq!(v3.2, 3.0);
    }

    #[test]
    fn vec_mul_f64() {
        let v1 = Tup(1.0, 2.0, 3.0);
        let a: f64 = 3.0;

        let v3 = v1 * a;
        assert_eq!(v3.0, 3.0);
        assert_eq!(v3.1, 6.0);
        assert_eq!(v3.2, 9.0);
    }

    #[test]
    fn vec_dot() {
        let v1 = Tup(1.0, 2.0, 3.0);
        let v2 = Tup(2.0, 3.0, 4.0);
        let a: f64 = v1.dot(v2);
        assert_eq!(a, 20.0);
    }

    #[test]
    fn vec_cross() {
        let v1 = Tup(1.0, 2.0, 3.0);
        let v2 = Tup(2.0, 3.0, 4.0);
        let a = v1.cross(v2);
        assert_eq!(a.0, -1.0);
        assert_eq!(a.1, 2.0);
        assert_eq!(a.2, -1.0);
    }

    #[test]
    fn vec_plus_equals() {
        let mut v1 = Tup(1.0, 2.0, 3.0);
        let v2 = Tup(2.0, 3.0, 4.0);
        v1 += v2;
        assert_eq!(v1.0, 3.0);
        assert_eq!(v1.1, 5.0);
        assert_eq!(v1.2, 7.0);
    }

    #[test]
    fn cartesian_to_boyer_lindquist() {
        let v = Tup(3., 4., 5.);
        let a = 0.999;
        let v_bl = v.cartesian_to_boyer_lindquist(a);
        // let expected_bl_coords = Tup(7.071, 1.107, 1.030);
        println!("v{:?}", v);
        println!("v_bl{:?}", v_bl);
        // println!("v_sp{:?}", v_bl.spherical_to_cartesian());
        println!("v{:?}", v_bl.boyer_lindquist_to_cartesian(a));
    }
}
